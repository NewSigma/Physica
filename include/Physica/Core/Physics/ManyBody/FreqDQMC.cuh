/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/UnitMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.cuh"

namespace Physica {
    template<Scalar T>
    class device_obj<FreqDQMC<T>> {
        using host_obj = FreqDQMC<T>;
        using This = device_obj<host_obj>;
        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;

        template<Scalar U>
        using HostMatrix = DenseMatrix<U, MatrixOption::Col, Dynamic, Dynamic, PageLockedAllocator<T>>;
    private:
        Array<device_obj<DenseLU<T, false>>, 2> lu;
        ActionMatrix<T> action;
        device_obj<MatrixND<T>> linearRHS;
        HostMatrix<T> detBuffer;
        device_obj<VectorND<T>> diag;

        Array<device_obj<MatrixND<Tr>>> greensD;
        Array<MatrixND<Tr>, 2> greensH;
        float64 lnAbsDet = Trv::nan();
        float64 sign = 1;
    public:
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
        void step(bool warmup = false);
        template<RNG R>
        void step_for(int numStep);

        [[nodiscard]] Vector2D<float64> promoteDet(const device_obj<MatrixND<T>>& matrixLU);
        void traceGreen(int spin);
        /* Getters */
        [[nodiscard]] int getNumSite() const noexcept { return action.getNumSite(); }
        [[nodiscard]] int getNumFreq() const noexcept { return action.getNumFreq(); }
        [[nodiscard]] const auto& getGreens() noexcept { return greensH; }
        [[nodiscard]] Trv getSign() const noexcept { return Trv(sign); }
        [[nodiscard]] Trv getRSign() const noexcept { return getSign(); }
    private:
        [[nodiscard]] Vector2D<float64> calcDet();
        void calcGreen();
    };

    template<Scalar T>
    device_obj<FreqDQMC<T>>::device_obj(const HubbardParams<Tr>& params, Trv freqDensity)
            : action(params, host_obj::calcFreqCutoff(params.getBeta(), freqDensity))
            , greensD(2, params.getNumSite())
            , greensH(2, params.getNumSite()) {
        size_t size = action.getOrder();
        for (auto& spinLU : lu)
            spinLU.resize(size);
        linearRHS.resize(size);
        detBuffer.resize(size);
        diag.resize(size);

        action.assign_kinetic(detBuffer);
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<FreqDQMC<T>>::step_random() {
        action.template randAuxField<R>();
        auto [lnAD, sgnD] = calcDet();
        lnAbsDet = lnAD;
        sign = sgnD;
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<FreqDQMC<T>>::step(bool warmup) {
        assert(lnAbsDet.isFinite() && "[Error]: Should random initialize before monte carlo step");
        const int site = std::uniform_int_distribution<int>(0, getNumSite() - 1)(R::getInstance());
        const int freq = std::uniform_int_distribution<int>(0, getNumFreq() * 2 - 1)(R::getInstance());
        const T save = action.template randAuxField<R>(freq, site);

        auto [lnAD, sgnD] = calcDet();
        const bool accept = float64::template random_uniform<R>() < exp(lnAD - lnAbsDet);
        if (accept) {
            lnAbsDet = lnAD;
            sign = sgnD;
            if (!warmup)
                calcGreen();
        }
        else
            action.getAuxField()(freq, site) = save;
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<FreqDQMC<T>>::step_for(int numStep) {
        for (int i = 0; i < numStep; ++i)
            step<R>(true);
        calcGreen();
    }
    /**
     * float64 is necessary for determinants
     */
    template<Scalar T>
    auto device_obj<FreqDQMC<T>>::promoteDet(const device_obj<MatrixND<T>>& matrixLU) -> Vector2D<float64> {
        auto kernel = [matrixLU_ = asStruct(matrixLU), diag_ = asStruct(diag)] __device__() mutable {
            const auto& matrixLU = matrixLU_.getDerived();
            auto& diag = diag_.getDerived();
            unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i < diag.getLength())
                diag[i] = matrixLU(i, i);
        };

        constexpr int WarpSize = CUDADevAttr::WarpSize;
        int numThread = WarpSize;
        int numBlock = (diag.getLength() + (WarpSize - 1)) / WarpSize;
        CUDAExecutor::launch(kernel, KernelConfig(numBlock, numThread));
        DiagMatrix<cfloat64> mat(diag.toHost());
        return {mat.lnAbsDet(), mat.sgndet().real()};
    }

    template<Scalar T>
    void device_obj<FreqDQMC<T>>::traceGreen(int spin) {
        auto kernel = [linearRHS_ = asStruct(linearRHS),
                       green = asStruct(greensD[spin]),
                       numSite = getNumSite(),
                       numFreq = getNumFreq(),
                       repBeta = reciprocal(action.getParams().getBeta())] __device__() mutable {
            const auto& linearRHS = linearRHS_.getDerived();
            unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;
            unsigned int col = blockIdx.y;
            Tr elem = 0;
            for (int _ = 0, offset = 0; _ < numFreq * 2; ++_) {
                elem += linearRHS(offset + row, offset + col).real();
                offset += numSite;
            }

            elem *= repBeta;
            if (row == col)
                elem += Tr(0.5);

            green.getDerived()(row, col) = elem;
        };

        int numSite = getNumSite();
        uint32_t numThread = std::min<uint32_t>(numSite, CUDADevAttr::WarpSize);
        uint32_t numBlockX = (numSite + numThread - 1) / numThread;
        uint32_t numBlockY = numSite;
        CUDAExecutor::launch(kernel, KernelConfig({numBlockX, numBlockY}, numThread));

        greensD[spin].toHostAsync(greensH[spin]);
    }

    template<Scalar T>
    auto device_obj<FreqDQMC<T>>::calcDet() -> Vector2D<float64> {
        action.assign_potential(detBuffer);
        lu[0].compute(detBuffer);

        action.flip();
        action.assign_potential(detBuffer);
        lu[1].compute(detBuffer);

        float64 lnAD = 0, sgnD = 1;
        for (const auto& spinLU : lu) {
            auto [x, y] = promoteDet(spinLU.getMatrixLU());
            lnAD += x;
            sgnD *= y;
        }
        return {lnAD, sgnD};
    }

    template<Scalar T>
    void device_obj<FreqDQMC<T>>::calcGreen() {
        for (int spin : {0, 1}) {
            auto& spinLU = lu[spin];
            linearRHS = device_obj<UnitMatrix<Tr>>(action.getOrder());
            spinLU.solve(linearRHS);
            traceGreen(spin);
        }
        CUDAExecutor::wait();
    }
}
