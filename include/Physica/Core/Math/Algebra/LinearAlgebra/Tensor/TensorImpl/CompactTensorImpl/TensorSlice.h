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

#include "../CompactTensor.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/CompactMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/StridedMatrix.h"

namespace Physica {
    namespace Internal {
        template<class X, int DimR, int DimC>
        __host__ __device__ consteval bool isTensorSliceCompact() noexcept {
            constexpr auto Strides = std::remove_cvref_t<X>::getStrideAtCompile();
            constexpr int Major = std::min(DimR, DimC);
            constexpr int Minor = std::max(DimR, DimC);
            return (Strides[Minor] == 1) && (Major + 1 == Minor);
        }

        template<class X, int DimR, int DimC>
        using TensorSliceBase = std::conditional_t<isTensorSliceCompact<X, DimR, DimC>(),
                                                   CompactMatrix<TensorSlice<X, DimR, DimC>>,
                                                   StridedMatrix<TensorSlice<X, DimR, DimC>>>;
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    class TensorSlice<X, DimR, DimC> : public Internal::TensorSliceBase<X, DimR, DimC> {
        using This = TensorSlice<X, DimR, DimC>;
        using Base = Internal::TensorSliceBase<X, DimR, DimC>;
        using IndexType = std::remove_cvref_t<X>::IndexType;

        static_assert(DimR != DimC, "[Error]: DimR and DimC must be different");
        static_assert(DimR < std::remove_cvref_t<X>::ndim());
        static_assert(DimC < std::remove_cvref_t<X>::ndim());
    private:
        decay_rvalue_t<X> tensor;
        IndexType index;
    public:
        TensorSlice(X&& tensor, IndexVar auto... indices);
        TensorSlice(const This&) = default;
        TensorSlice(This&&) noexcept = default;
        ~TensorSlice() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(size_t row, size_t col);
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return tensor.dim(DimR); }
        [[nodiscard]] size_t getCol() const noexcept { return tensor.dim(DimC); }
        [[nodiscard]] size_t getOrder() const noexcept;
        [[nodiscard]] constexpr size_t getRowStride() const noexcept;
        [[nodiscard]] constexpr size_t getColStride() const noexcept;
        [[nodiscard]] auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept { return DimR < DimC ? MatrixMajor::Row : MatrixMajor::Col; }
        [[nodiscard]] __host__ __device__ consteval static size_t getRowStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept;
    };

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    TensorSlice<X, DimR, DimC>::TensorSlice(X&& tensor, IndexVar auto... indices) : tensor(std::forward<X>(tensor)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(indices < tensor.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    void TensorSlice<X, DimR, DimC>::resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow() && col == getCol());
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    size_t TensorSlice<X, DimR, DimC>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return getRow();
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    constexpr size_t TensorSlice<X, DimR, DimC>::getRowStride() const noexcept {
        return tensor.getStride(DimR);
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    constexpr size_t TensorSlice<X, DimR, DimC>::getColStride() const noexcept {
        return tensor.getStride(DimC);
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    auto TensorSlice<X, DimR, DimC>::data_handle(this auto&& self) noexcept {
        auto idx = self.index;
        idx[DimR] = 0;
        idx[DimC] = 0;
        return self.tensor.data_ptr(idx);
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ consteval size_t TensorSlice<X, DimR, DimC>::getRowStrideAtCompile() noexcept {
        return std::remove_cvref_t<X>::getStrideAtCompile()[DimR];
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ consteval size_t TensorSlice<X, DimR, DimC>::getColStrideAtCompile() noexcept {
        return std::remove_cvref_t<X>::getStrideAtCompile()[DimC];
    }
}
