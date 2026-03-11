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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.cuh"
#include "DQMCImpl/ActionMatrix.cuh"

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
        using HostMatrix = DenseMatrix<U, MatrixMajor::Col, Dynamic, Dynamic, PageLockedAllocator<U>>;

        static_assert(T::Prec == Float32, "[Info]: Suggest FP32");
    private:
        Array<device_obj<DenseLU<T, false>>, 2> lu;
        device_obj<ActionMatrix<T>> action;
        device_obj<MatrixND<T>> forceBuffer;
        device_obj<MatrixND<T>> solBuffer;
        Trv correction;

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
        [[nodiscard]] VectorND<Trv> makeDefaultMass() const;

        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Vector auto& pos);
        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Cell& cell);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Vector auto& pos, Vector auto& result);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Cell& cell, Vector auto& result);
        /* Getters */
        [[nodiscard]] auto&& getAuxField(this auto&& self) noexcept { return self.action.getAuxField(); }
        [[nodiscard]] int getNumSite() const noexcept { return action.getNumSite(); }
        [[nodiscard]] int getNumFreq() const noexcept { return action.getNumFreq(); }
        [[nodiscard]] int getMaxBoson() const noexcept { return action.getMaxBoson(); }
        [[nodiscard]] const auto& getGreens() noexcept { return greensH; }
        [[nodiscard]] Trv getSign() const noexcept { return Trv(sign); }
        [[nodiscard]] Trv getRSign() const noexcept { return getSign(); }
        [[nodiscard]] Trv getRepelU() const noexcept { return action.getRepelU(); }
        [[nodiscard]] Trv getBeta() const noexcept { return action.getBeta(); }
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
            : action(params, ElasticDQMC<Trv>::calcFreqCutoff(params.getBeta(), freqDensity), maxBoson)
            , correction(ElasticDQMC<Trv>::calcLocalCorrection(params.getBeta(), params.getRepelU(), params.getChemMu(), getNumFreq()))
            , greensD(2, params.getNumSite())
            , greensH(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
        forceBuffer.resize(getAuxField());
        solBuffer.resize(size);
    }

    template<Scalar T>
    device_obj<FreqDQMC<T>>::device_obj(const HubbardParams<Tr>& params, Trv freqDensity)
            : action(params, ElasticDQMC<Trv>::calcFreqCutoff(params.getBeta(), freqDensity))
            , correction(ElasticDQMC<Trv>::calcLocalCorrection(params.getBeta(), params.getRepelU(), params.getChemMu(), getNumFreq()))
            , greensD(2, params.getNumSite())
            , greensH(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
        forceBuffer.resize(getAuxField());
        solBuffer.resize(size);
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
        check(cudaMemcpyAsync(init.data(), getAuxField().data(), sizeof(Tr) * init.getLength(), cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
        CUDAExecutor::wait();
        hmc.setInitPosition(std::move(init));
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    auto device_obj<FreqDQMC<T>>::step(HamiltonMC<Tr>& hmc) -> Trv {
        Trv acceptR = hmc.template step<R, P>(*this);
        check(cudaMemcpyAsync(getAuxField().data(), hmc.getSample().data(), sizeof(Tr) * hmc.getDOF(), cudaMemcpyHostToDevice, CUDAContext::getInstance()));

        auto [lnAD, sgnD] = calcDet();
        lnWeight = lnAD;
        sign = sgnD;
        calcGreen();
        return acceptR;
    }

    template<Scalar T>
    auto device_obj<FreqDQMC<T>>::makeDefaultMass() const -> VectorND<Trv> {
        VectorND<Trv> result(getAuxField().getSize() * 2, 1);
        result.tail(2) *= Trv(2);
        return result;
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto device_obj<FreqDQMC<T>>::potentialV(const Vector auto& pos) -> Trv {
        assert(pos.getLength() == getAuxField().getSize() * 2 && "[Error]: Real matrix contains 2x number of elements of complex matrix");
        check(cudaMemcpyAsync(getAuxField().data(), pos.data(), sizeof(Tr) * pos.getLength(), cudaMemcpyHostToDevice, CUDAContext::getInstance()));

        MatrixND<T> buffer;
        buffer.resize(getAuxField());
        buffer.read(pos);

        MatrixND<Tr> sqnorms = buffer.squaredNorms();
        sqnorms *= reciprocal(getBetaU());
        sqnorms.row(0) *= Trv(0.5);
        return sqnorms.sum() - calcDet()[0];
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
        check(cudaMemcpyAsync(getAuxField().data(), pos.data(), sizeof(Tr) * pos.getLength(), cudaMemcpyHostToDevice, CUDAContext::getInstance()));
        forceBuffer.zeros();

        for (int spin : {0, 1}) {
            auto& spinLU = lu[spin];
            spinLU.compute(action);
            solBuffer = spinLU.inv();

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
            int numThread = std::min<int>(numSite, CUDADevAttr::DefaultThreadsPerBlock);
            int numBlockX = getMaxBoson();
            CUDAExecutor::launch(collectForce, KernelConfig(numBlockX, numThread));

            action.flip();
        }

        MatrixND<T> force;
        force.resize(getAuxField());
        force.read(pos);
        force *= Trv(-2) / getBetaU();
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
        lu[0].compute(action);

        action.flip();
        lu[1].compute(action);

        Trv lnAD = 0, sgnD = 1;
        for (auto& spinLU : lu) {
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
            solBuffer = lu[spin].inv();
            traceGreen(spin);
        }
        CUDAExecutor::wait();
    }

    template<Scalar T>
    void device_obj<FreqDQMC<T>>::traceGreen(int spin) {
        auto kernel = [solBuffer_ = asStruct(solBuffer),
                       green = asStruct(greensD[spin]),
                       numSite = getNumSite(),
                       size = 2 * getNumFreq(),
                       correction = correction] __device__() mutable {
            const auto& solBuffer = solBuffer_.getDerived();
            unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
            unsigned int col = blockIdx.y;
            Tr elem = 0;
            for (int _ = 0, offset = 0; _ < size; ++_) {
                elem += solBuffer[offset + row, offset + col].real();
                offset += numSite;
            }

            if (row == col)
                elem += Tr(0.5) + correction;

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
        return getBeta() * getRepelU();
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<device_obj<FreqDQMC<T>>> : public Traits<FreqDQMC<T>> {};
}
