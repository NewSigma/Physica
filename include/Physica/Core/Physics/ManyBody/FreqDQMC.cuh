/*
 * Copyright 2025-2026 Weibo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#include "FreqDQMC.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/IdentityMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.cuh"

namespace Physica {
    template<Scalar T>
    class device_obj<FreqDQMC<T>> {
        using host_obj = FreqDQMC<T>;
        using This = device_obj<host_obj>;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        using Cell = MDCell<Tr, 1>;
        template<Scalar U>
        using HostMatrix = DenseMatrix<U, MatrixOption::Col, Dynamic, Dynamic, PageLockedAllocator<U>>;
    private:
        Array<device_obj<DenseLU<T, false>>, 2> lu;
        device_obj<MatrixND<T>> forceBuffer{};
        device_obj<MatrixND<T>> solBuffer;
        HostMatrix<T> detBuffer;
        ActionMatrix<T> action;

        Array<device_obj<MatrixND<Tr>>, 2> greensD;
        Array<MatrixND<Tr>, 2> greensH;
        Trv lnWeight = Trv::nan();
        Trv sign = 1;
    public:
        device_obj(const HubbardParams<Tr>& params, Trv freqDensity, int maxBoson);
        device_obj(const HubbardParams<Tr>& params, Trv freqDensity);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<RNG R>
        void step_random();
        template<RNG R>
        bool step(bool warmup = false);
        template<RNG R>
        void step_for(int numStep);

        template<RNG R>
        void step_random(HamiltonMC<Tr>& hmc);
        template<RNG R, ExecutePolicy P = Sequential>
        Trv step(HamiltonMC<Tr>& hmc);

        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Vector auto& pos);
        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Cell& cell);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Vector auto& pos, Vector auto& result);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Cell& cell, Vector auto& result);
        /* Getters */
        [[nodiscard]] const auto& getParams() const noexcept { return action.getParams(); }
        [[nodiscard]] auto& getAuxField() noexcept { return action.getAuxField(); }
        [[nodiscard]] int getNumSite() const noexcept { return action.getNumSite(); }
        [[nodiscard]] int getNumFreq() const noexcept { return action.getNumFreq(); }
        [[nodiscard]] int getMaxBoson() const noexcept { return action.getMaxBoson(); }
        [[nodiscard]] const auto& getGreens() noexcept { return greensH; }
        [[nodiscard]] Trv getSign() const noexcept { return Trv(sign); }
        [[nodiscard]] Trv getRSign() const noexcept { return getSign(); }
    private:
        [[nodiscard]] Vector2D<Trv> calcDet();
        [[nodiscard]] Vector2D<Trv> calcLnWeight();
        void calcGreen();
        void traceGreen(int spin);
        /* Getters */
        [[nodiscard]] Trv getBetaU() const noexcept;
    };

    template<Scalar T>
    device_obj<FreqDQMC<T>>::device_obj(const HubbardParams<Tr>& params, Trv freqDensity, int maxBoson)
            : action(params, host_obj::calcFreqCutoff(params.getBeta(), freqDensity), maxBoson)
            , greensD(2, params.getNumSite())
            , greensH(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
        forceBuffer.resize(getAuxField());
        solBuffer.resize(size);
        detBuffer.resize(size);

        action.assign_kinetic(detBuffer);
    }

    template<Scalar T>
    device_obj<FreqDQMC<T>>::device_obj(const HubbardParams<Tr>& params, Trv freqDensity)
            : action(params, host_obj::calcFreqCutoff(params.getBeta(), freqDensity))
            , greensD(2, params.getNumSite())
            , greensH(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
        forceBuffer.resize(getAuxField());
        solBuffer.resize(size);
        detBuffer.resize(size);

        action.assign_kinetic(detBuffer);
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<FreqDQMC<T>>::step_random() {
        action.template randAuxField<R>();
        auto [lnW, sgnD] = calcLnWeight();
        lnWeight = lnW;
        sign = sgnD;
    }

    template<Scalar T>
    template<RNG R>
    bool device_obj<FreqDQMC<T>>::step(bool warmup) {
        assert(lnWeight.isFinite() && "[Error]: Should random initialize before monte carlo step");
        const int site = std::uniform_int_distribution<int>(0, getNumSite() - 1)(R::getInstance());
        const int freq = std::uniform_int_distribution<int>(0, getMaxBoson() - 1)(R::getInstance());
        const T save = action.template randAuxField<R>(freq, site);

        auto [lnW, sgnD] = calcLnWeight();
        const bool accept = Trv::template random_uniform<R>() < exp(lnW - lnWeight);
        if (accept) {
            lnWeight = lnW;
            sign = sgnD;
            if (!warmup)
                calcGreen();
        }
        else
            getAuxField()[freq, site] = save;
        return accept;
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<FreqDQMC<T>>::step_for(int numStep) {
        for (int i = 0; i < numStep; ++i)
            step<R>(true);
        calcGreen();
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<FreqDQMC<T>>::step_random(HamiltonMC<Tr>& hmc) {
        action.template randAuxField<R>();

        VectorND<Tr> init(getAuxField().getSize() * 2);
        init.read(getAuxField());
        hmc.setInitPosition(std::move(init));
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    auto device_obj<FreqDQMC<T>>::step(HamiltonMC<Tr>& hmc) -> Trv {
        Trv acceptR = hmc.template step<R, P>(*this);
        getAuxField().read(hmc.getSample());

        auto [lnAD, sgnD] = calcDet();
        lnWeight = lnAD;
        sign = sgnD;
        calcGreen();
        return acceptR;
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto device_obj<FreqDQMC<T>>::potentialV(const Vector auto& pos) -> Trv {
        assert(pos.getLength() == getAuxField().getSize() * 2 && "[Error]: Real matrix contains 2x number of elements of complex matrix");
        getAuxField().read(pos);

        MatrixND<Tr> buffer = getAuxField().squaredNorms();
        buffer *= reciprocal(getBetaU());
        buffer.row(0) *= Trv(0.5);
        return buffer.sum() - calcDet()[0];
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto device_obj<FreqDQMC<T>>::potentialV(const Cell& cell) -> Trv {
        return potentialV<P>(cell.getPos().flatten());
    }

    template<Scalar T>
    template<ExecutePolicy P>
    void device_obj<FreqDQMC<T>>::forceAsync(const Vector auto& pos, Vector auto& result) {
        assert(result.getLength() == pos.getLength());
        assert(result.getLength() == getAuxField().getSize() * 2);
        getAuxField().read(pos);
        forceBuffer.zeros();

        for (int spin : {0, 1}) {
            auto& spinLU = lu[spin];
            action.assign_potential(detBuffer);
            spinLU.compute(detBuffer);

            solBuffer = device_obj<IdentityMatrix<Tr>>(action.getOrder());
            spinLU.solve(solBuffer);

            Trv factor = spin == 0 ? 1.0 : -1.0;
            int numSite = getNumSite();
            auto collectForce = [spinF_ = asStruct(forceBuffer), inv_ = asStruct(solBuffer), factor, numSite] __device__() mutable {
                const auto& inv = inv_.getDerived();
                auto& spinF = spinF_.getDerived();
                unsigned int delta = blockIdx.x;
                unsigned int size = gridDim.x;
                ThreadBlock block{};

                bool elastic = delta == 0;
                if (elastic) {
                    for (unsigned int i = 0; i < size; ++i) {
                        const auto diag = inv.transpose().block(i * numSite, numSite, i * numSite, numSite).diag();
                        (diag.reals() * factor).assign_add(spinF.row(0), block);
                    }
                }
                else {
                    unsigned int r = 0;
                    unsigned int c = delta;
                    auto row = spinF.row(delta);
                    for (; c < size; ++r, ++c) {
                        const auto upper = inv.transpose().block(r * numSite, numSite, c * numSite, numSite).diag();
                        const auto lower = inv.transpose().block(c * numSite, numSite, r * numSite, numSite).diag();
                        (upper * factor).assign_add(row, block);
                        (lower.conjugate() * factor).assign_add(row, block);
                    }
                }
            };
            int numThread = std::min(numSite, device_obj<VectorND<T>>::MaxThreadsPerBlock);
            int numBlockX = getMaxBoson();
            CUDAExecutor::launch(collectForce, KernelConfig(numBlockX, numThread));

            action.flip();
        }

        MatrixND<T> force = -getAuxField() * (Trv(2) / getBetaU());
        force.row(0) *= Trv(0.5);
        force += forceBuffer.toHost();
        result.read(force);
    }

    template<Scalar T>
    template<ExecutePolicy P>
    void device_obj<FreqDQMC<T>>::forceAsync(const Cell& cell, Vector auto& result) {
        forceAsync<P>(cell.getPos().flatten(), result);
    }

    template<Scalar T>
    auto device_obj<FreqDQMC<T>>::calcDet() -> Vector2D<Trv> {
        action.assign_potential(detBuffer);
        lu[0].compute(detBuffer);

        action.flip();
        action.assign_potential(detBuffer);
        lu[1].compute(detBuffer);

        Trv lnAD = 0, sgnD = 1;
        for (const auto& spinLU : lu) {
            lnAD += spinLU.lnAbsDet();
            sgnD *= spinLU.sgndet().real();
        }
        return {lnAD, sgnD};
    }

    template<Scalar T>
    auto device_obj<FreqDQMC<T>>::calcLnWeight() -> Vector2D<Trv> {
        Trv betaU = getBetaU();
        auto [lnAD, sgnD] = calcDet();
        lnAD = lnAD - ln1pexp(lncosh(getAuxField().row(0).reals()) + fma(betaU, Trv(-0.5), MathConst<Trv>::ln2)).sum();
        if (getMaxBoson() > 1) {
            lnAD -= ln1pexp(lncosh(getAuxField().bottomRows(1).flatten().reals()) + fma(betaU, Trv(-0.25), MathConst<Trv>::ln2)).sum()
                  + ln1pexp(lncosh(getAuxField().bottomRows(1).flatten().imags()) + fma(betaU, Trv(-0.25), MathConst<Trv>::ln2)).sum();
        }
        return {lnAD, sgnD};
    }

    template<Scalar T>
    void device_obj<FreqDQMC<T>>::calcGreen() {
        for (int spin : {0, 1}) {
            auto& spinLU = lu[spin];
            solBuffer = device_obj<IdentityMatrix<Tr>>(action.getOrder());
            spinLU.solve(solBuffer);
            traceGreen(spin);
        }
        CUDAExecutor::wait();
    }

    template<Scalar T>
    void device_obj<FreqDQMC<T>>::traceGreen(int spin) {
        auto kernel = [solBuffer_ = asStruct(solBuffer),
                       green = asStruct(greensD[spin]),
                       numSite = getNumSite(),
                       size = 2 * getNumFreq()] __device__() mutable {
            const auto& solBuffer = solBuffer_.getDerived();
            unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
            unsigned int col = blockIdx.y;
            Tr elem = 0;
            for (int _ = 0, offset = 0; _ < size; ++_) {
                elem += solBuffer[offset + row, offset + col].real();
                offset += numSite;
            }

            if (row == col)
                elem += Tr(0.5);

            green.getDerived()[row, col] = elem;
        };

        int numSite = getNumSite();
        uint32_t numThread = std::min<uint32_t>(numSite, CUDADevAttr::WarpSize);
        uint32_t numBlockX = (numSite + numThread - 1) / numThread;
        uint32_t numBlockY = numSite;
        CUDAExecutor::launch(kernel, KernelConfig({numBlockX, numBlockY}, numThread));

        greensD[spin].toHostAsync(greensH[spin]);
    }

    template<Scalar T>
    auto device_obj<FreqDQMC<T>>::getBetaU() const noexcept -> Trv {
        return getParams().getBeta() * getParams().getRepelU();
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<device_obj<FreqDQMC<T>>> : public Traits<FreqDQMC<T>> {};
}
