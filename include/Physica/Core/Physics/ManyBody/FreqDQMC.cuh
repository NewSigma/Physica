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
        HubbardParams<Tr> params;
        Array<device_obj<DenseLU<T, false>>, 2> lu;
        device_obj<ActionMatrix<T>> action;
        device_obj<MatrixND<T>> forceBuffer;
        device_obj<MatrixND<T>> solBuffer;
        HamiltonMC<Tr> hmc;
        Trv correction;

        Array<device_obj<MatrixND<Tr>>, 2> greensD;
        Array<MatrixND<Tr>, 2> greensH;
        Trv lnWeight = Trv::nan();
        Trv sign = 1;
    public:
        device_obj(HubbardParams<Tr> params_, Trv freqDensity, int maxBoson);
        device_obj(HubbardParams<Tr> params_, Trv freqDensity);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        template<RNG R>
        void step_random();
        template<RNG R, ExecutePolicy P = Sequential>
        Trv step();
        template<RNG R, ExecutePolicy P = Sequential>
        void step_for(int numStep);

        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Vector auto& pos);
        template<ExecutePolicy P = Sequential>
        [[nodiscard]] Trv potentialV(const Cell& cell);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Vector auto& pos, Vector auto& result);
        template<ExecutePolicy P = Sequential>
        void forceAsync(const Cell& cell, Vector auto& result);

        Tr calcBerry();
        /* Getters */
        [[nodiscard]] auto&& getAuxField(this auto&& self) noexcept { return self.action.getAuxField(); }
        [[nodiscard]] int getNumSite() const noexcept { return action.getNumSite(); }
        [[nodiscard]] int getNumFreq() const noexcept { return action.getNumFreq(); }
        [[nodiscard]] int getMaxBoson() const noexcept { return action.getMaxBoson(); }
        [[nodiscard]] auto&& getHMC(this auto&& self) noexcept { return self.hmc; }
        [[nodiscard]] const auto& getGreens() noexcept { return greensH; }
        [[nodiscard]] Trv getSign() const noexcept { return Trv(sign); }
        [[nodiscard]] Trv getRSign() const noexcept { return getSign(); }
        [[nodiscard]] Trv getRepelU() const noexcept { return action.getRepelU(); }
        [[nodiscard]] Trv getBeta() const noexcept { return action.getBeta(); }
        /* Static members */
        [[nodiscard]] consteval static bool isPeriodBoundary() noexcept { return host_obj::isPeriodBoundary(); }
    private:
        [[nodiscard]] Vector2D<Trv> calcDet();
        void calcGreen();
        /* Getters */
        [[nodiscard]] Trv getBetaU() const noexcept;
    };

    template<Scalar T>
    device_obj<FreqDQMC<T>>::device_obj(HubbardParams<Tr> params_, Trv freqDensity, int maxBoson)
            : params(std::move(params_))
            , action(params, ElasticDQMC<Trv>::calcFreqCutoff(params.getBeta(), freqDensity), maxBoson)
            , hmc(host_obj::makeDefaultMass(getAuxField()))
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
    device_obj<FreqDQMC<T>>::device_obj(HubbardParams<Tr> params_, Trv freqDensity)
            : params(std::move(params_))
            , action(params, ElasticDQMC<Trv>::calcFreqCutoff(params.getBeta(), freqDensity))
            , hmc(host_obj::makeDefaultMass(getAuxField()))
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
        action.template random_normal<R>();

        VectorND<Tr> init(getAuxField().getSize() * 2);
        check(cudaMemcpyAsync(init.data(), getAuxField().data(), sizeof(Tr) * init.getLength(), cudaMemcpyDeviceToHost, CUDAContext::getInstance()));
        CUDAExecutor::wait();
        hmc.setInitPosition(std::move(init));
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    auto device_obj<FreqDQMC<T>>::step() -> Trv {
        Trv acceptR = hmc.template step<R, P>(*this);
        bool hasAuxField = !getBetaU().isSubNormal();
        if (hasAuxField)
            check(cudaMemcpyAsync(getAuxField().data(), hmc.getSample().data(), sizeof(Tr) * hmc.getDOF(), cudaMemcpyHostToDevice, CUDAContext::getInstance()));
        else
            getAuxField().zeros();

        auto [lnAD, sgnD] = calcDet();
        lnWeight = lnAD;
        sign = sgnD;
        calcGreen();
        return acceptR;
    }

    template<Scalar T>
    template<RNG R, ExecutePolicy P>
    void device_obj<FreqDQMC<T>>::step_for(int numStep) {
        for (int _ = 0; _ < numStep; ++_)
            step<R, P>();
    }

    template<Scalar T>
    template<ExecutePolicy P>
    auto device_obj<FreqDQMC<T>>::potentialV(const Vector auto& pos) -> Trv {
        assert(pos.getLength() == getAuxField().getSize() * 2 && "[Error]: Real matrix contains 2x number of elements of complex matrix");
        check(cudaMemcpyAsync(getAuxField().data(), pos.data(), sizeof(Tr) * pos.getLength(), cudaMemcpyHostToDevice, CUDAContext::getInstance()));
        bool noAuxField = getBetaU().isSubNormal();
        if (noAuxField)
            return -calcDet()[0];

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

            int numSite = getNumSite();
            int numThread = std::min<int>(numSite, CUDADevAttr::DefaultThreadsPerBlock);
            int numBlockX = getMaxBoson();
            CUDAExecutor::launch([spinF_ = asStruct(forceBuffer),
                                 inv_ = asStruct(solBuffer),
                                 factor = Trv(spin == 0 ? 1 : -1),
                                 numSite,
                                 numFreq2 = getNumFreq() * 2] __device__() mutable {
                const auto invT = inv_.getDerived().transpose();
                auto& spinF = spinF_.getDerived();
                const int freq = int(blockIdx.x);
                for (int site = int(threadIdx.x); site < numSite; site += int(blockDim.x)) {
                    const auto block = invT.block(site * numFreq2, numFreq2, site * numFreq2, numFreq2);
                    if (freq == 0)
                        spinF[0, site] += block.diag().sum().real() * factor;
                    else
                        spinF[freq, site] += (block.diag(freq).sum() + block.diag(-freq).conjugate().sum()) * factor;
                }
            }, KernelConfig(numBlockX, numThread));
            action.flip();
        }

        MatrixND<T> force(getAuxField().getRow(), getAuxField().getCol());
        bool noAuxField = getBetaU().isSubNormal();
        if (noAuxField)
            force.zeros();
        else {
            force.read(pos);
            force *= Trv(-2) / getBetaU();
            force.row(0) *= Trv(0.5);
        }
        force += forceBuffer.toHost();
        result.read(force);
    }

    template<Scalar T>
    template<ExecutePolicy P>
    void device_obj<FreqDQMC<T>>::forceAsync(const Cell& cell, Vector auto& result) {
        forceAsync<P>(cell.getPos().flatten(), result);
    }

    template<Scalar T>
    auto device_obj<FreqDQMC<T>>::calcBerry() -> Tr {
        const int numSite = getNumSite();
        const int numFreq2 = 2 * getNumFreq();
        device_obj<VectorND<Tr>> ptrace(numFreq2);
        Tr result = 0;
        for (auto& spinLU : lu) {
            solBuffer = spinLU.inv();

            ptrace.zeros();
            auto kernel = [inv_ = asStruct(solBuffer), ptrace_ = asStruct(ptrace), numSite] __device__() mutable {
                const auto& inv = inv_.getDerived();
                auto& ptrace = ptrace_.getDerived();

                int numFreq2 = ptrace.getLength();
                int i = blockIdx.x * blockDim.x + threadIdx.x;
                if (i < numFreq2) {
                    for (int site = 0; site < numSite; ++site)
                        ptrace[i] += inv.block(site * numFreq2, numFreq2, site * numFreq2, numFreq2).diag().imags().calc(i);
                }
            };
            CUDAExecutor::launch(kernel, ptrace.makeKernelConfig());
            result -= ptrace * action.getMatsubara().diag();
        }
        const Tr background = reciprocal(Trv(1) + square(divide(params.calcBetaMu(), action.getMatsubara().diag().head(getNumFreq())))).sum() * (4 * numSite);
        return result - background;
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
    void device_obj<FreqDQMC<T>>::calcGreen() {
        for (int spin : {0, 1}) {
            solBuffer = lu[spin].inv();
            auto kernel = [solBuffer_ = asStruct(solBuffer),
                        green_ = asStruct(greensD[spin]),
                        numSite = getNumSite(),
                        size = 2 * getNumFreq(),
                        correction = correction] __device__() mutable {
                const auto& solBuffer = solBuffer_.getDerived();
                auto& green = green_.getDerived();
                unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
                unsigned int col = blockIdx.y;
                if (row >= green.getRow())
                    return;

                Tr elem = 0;
                for (int i = 0, offset = 0; i < size; ++i) {
                    elem += solBuffer[offset + row, offset + col].real();
                    offset += numSite;
                }

                if (row == col)
                    elem += Tr(0.5) + correction;

                green[row, col] = elem;
            };

            int numSite = getNumSite();
            uint32_t numThread = std::min<uint32_t>(numSite, CUDADevAttr::WarpSize);
            uint32_t numBlockX = (numSite + numThread - 1) / numThread;
            uint32_t numBlockY = numSite;
            CUDAExecutor::launch(kernel, KernelConfig({numBlockX, numBlockY}, numThread));

            greensD[spin].toHostAsync(greensH[spin]);
        }
        CUDAExecutor::wait();
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
