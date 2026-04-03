/*
 * Copyright 2026 Weibo He.
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

#include "ActionMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiagMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/IdentityMatrix.cuh"

namespace Physica {
    template<Scalar T>
    class device_obj<ActionMatrix<T>> : public device_obj<RValueMatrix<ActionMatrix<T>>> {
        using host_obj = ActionMatrix<T>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;

        using Tr = T::RealType;
        using Tv = T::ValueType;
        using Trv = Tr::ValueType;
        static_assert(T::isComplex(), "[Error]: Action is complex");
    private:
        device_obj<DiagMatrix<Tr>> matsubara;
        device_obj<MatrixND<T>> auxField;
        device_obj<MatrixND<Tr>> hoppingT;
        Tr repelU;
        Tr chemMu;
        Tr beta;
    public:
        device_obj(const HubbardParams<Tr>& params, int numFreq, int maxBoson);
        device_obj(const HubbardParams<Tr>& params, int numFreq);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto&& target) const;
        void assign_kinetic(Matrix auto&& target) const;
        void assign_potential(Matrix auto&& target) const;

        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;

        [[nodiscard]] __host__ __device__ auto&& transpose(this auto&&) noexcept;

        void flip();
        template<RNG R>
        void randAuxField();
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return matsubara.getOrder() * getNumSite(); }
        [[nodiscard]] size_t getRow() const noexcept { return getOrder(); }
        [[nodiscard]] size_t getCol() const noexcept { return getOrder(); }
        [[nodiscard]] auto&& getAuxField(this auto&&) noexcept;
        [[nodiscard]] int getNumFreq() const noexcept { return matsubara.getOrder() / 2; }
        [[nodiscard]] int getNumSite() const noexcept { return auxField.getCol(); }
        [[nodiscard]] int getMaxBoson() const noexcept { return auxField.getRow(); }
        [[nodiscard]] Tr getRepelU() const noexcept { return repelU; }
        [[nodiscard]] Tr getBeta() const noexcept { return beta; }
    };

    template<Scalar T>
    device_obj<ActionMatrix<T>>::device_obj(const HubbardParams<Tr>& params, int numFreq, int maxBoson)
            : matsubara(numFreq * 2)
            , auxField(maxBoson, params.getNumSite())
            , hoppingT(params.getHoppingMatrix())
            , repelU(params.getRepelU())
            , chemMu(params.getChemMu())
            , beta(params.getBeta()) {
        assert(numFreq > 0);
        assert((1 <= maxBoson) && (maxBoson <= 2 * numFreq) && "[Error]: maxBoson out of range");
        VectorND<Tr> diag(matsubara.getOrder());
        for (int k = 0; k < diag.getLength(); ++k) {
            int m = k - numFreq;
            diag[k] = Trv(2 * m + 1);
        }
        diag *= MathConst<Trv>::pi;
        diag.toDevice(matsubara.diag());
    }

    template<Scalar T>
    device_obj<ActionMatrix<T>>::device_obj(const HubbardParams<Tr>& params, int numFreq) : device_obj(params, numFreq, numFreq * 2) {}

    template<Scalar T>
    void device_obj<ActionMatrix<T>>::assign(Matrix auto&& target) const {
        // Note that we seldom change the kinetic part, separate them to customize the potential part.
        assign_kinetic(target);
        assign_potential(target);
    }

    template<Scalar T>
    void device_obj<ActionMatrix<T>>::assign_kinetic(Matrix auto&& target) const {
        kronecker(device_obj<IdentityMatrix<Trv>>(getNumSite()) * T(0, 1), matsubara).assign(target);
        kronecker(hoppingT * beta, device_obj<IdentityMatrix<Trv>>(matsubara.getOrder())).assign_add(target);
    }

    template<Scalar T>
    void device_obj<ActionMatrix<T>>::assign_potential(Matrix auto&& target) const {
        const unsigned int numThread = std::min<unsigned int>(getNumSite(), CUDADevAttr::DefaultThreadsPerBlock);
        const unsigned int numBlockX = 1;
        const unsigned int numBlockY = getNumFreq() * 2;
        const unsigned int numBlockZ = getNumFreq() * 2;
        CUDAExecutor::launch([auxField_ = asStruct(auxField),
                              target_ = asStruct(target),
                              maxBoson = getMaxBoson(),
                              shift = beta * fma(repelU, Tr(-0.5), chemMu)] __device__() mutable {
            int r = blockIdx.y;
            int c = blockIdx.z;
            int delta = std::abs(r - c);
            if (delta >= maxBoson)
                return;

            int site = threadIdx.x;
            int numFreq2 = gridDim.y;
            
            int offset = site * numFreq2;
            auto block = target_.getDerived().block(offset, numFreq2, offset, numFreq2);
            const auto& auxField = auxField_.getDerived();
            if (r > c)
                block[r, c] = auxField[delta, site];
            else if (r < c)
                block[r, c] = auxField[delta, site].conjugate();
            else
                block[r, r].real() = auxField[0, site].real() - shift;
        }, KernelConfig({numBlockX, numBlockY, numBlockZ}, numThread));
    }

    template<Scalar T>
    __device__ T device_obj<ActionMatrix<T>>::calc(size_t row, size_t col) const {
        assert(row < getRow());
        assert(col < getCol());
        const int numSite = getNumSite();
        const size_t rowSite = row / numSite;
        const size_t rowFreq = row % numSite;
        const size_t colSite = col / numSite;
        const size_t colFreq = col % numSite;
        bool diagFreq = rowFreq == colFreq;
        bool diagSite = rowSite == colSite;
        if (diagFreq) {
            Tr re = 0, im = 0;
            re = hoppingT[rowSite, colSite] * beta;
            if (diagSite) {
                const Tr shift = beta * fma(repelU, Tr(-0.5), chemMu);
                re += auxField[0, rowSite].real() - shift;
                im = matsubara.diag()[rowFreq];
            }
            return T(re, im);
        }

        if (diagSite) {
            bool upper = rowFreq > colFreq;
            auto delta = upper ? (rowFreq - colFreq) : (colFreq - rowFreq);
            if (delta < getMaxBoson()) {
                T aux = auxField[delta, rowSite];
                return upper ? aux : aux.conjugate();
            }
        }
        return 0;
    }

    template<Scalar T>
    __host__ __device__ auto&& device_obj<ActionMatrix<T>>::transpose(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self);
    }

    template<Scalar T>
    void device_obj<ActionMatrix<T>>::flip() {
        auxField = -auxField;
    }

    template<Scalar T>
    template<RNG R>
    void device_obj<ActionMatrix<T>>::randAuxField() {
        Tr factor = sqrt(repelU * beta);
        MatrixND<T> buffer;
        buffer.resize(auxField);
        buffer.template random_normal<R>();
        buffer.row(0) *= factor;
        if (getMaxBoson() > 1)
            buffer.bottomRows(1) *= factor / sqrt(Trv(2));
        buffer.toDevice(auxField);
    }

    template<Scalar T>
    auto&& device_obj<ActionMatrix<T>>::getAuxField(this auto&& self) noexcept {
        return std::forward<decltype(self)>(self).auxField;
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<device_obj<ActionMatrix<T>>> : public Traits<ActionMatrix<T>> {};
}
