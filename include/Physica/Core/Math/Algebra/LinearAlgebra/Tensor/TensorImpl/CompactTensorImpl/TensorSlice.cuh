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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/CompactMatrix.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/StridedMatrix.cuh"
#include "../CompactTensor.cuh"
#include "TensorSlice.h"

namespace Physica {
    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    class device_obj<TensorSlice<X, DimR, DimC>> : public device_obj<Internal::TensorSliceBase<X, DimR, DimC>> {
        using host_obj = TensorSlice<X, DimR, DimC>;
        using This = device_obj<host_obj>;
        using Base = device_obj<Internal::TensorSliceBase<X, DimR, DimC>>;
        using Ref = add_device_obj_t<X>;
        using IndexType = std::remove_cvref_t<X>::IndexType;

        static_assert(DimR != DimC, "[Error]: DimR and DimC must be different");
        static_assert(DimR < std::remove_cvref_t<X>::ndim());
        static_assert(DimC < std::remove_cvref_t<X>::ndim());
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<X>>> tensor;
        IndexType index;
    public:
        __host__ __device__ device_obj(Ref tensor, IndexVar auto... indices);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize(size_t row, size_t col);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return tensor.getDerived().dim(DimR); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return tensor.getDerived().dim(DimC); }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getRowStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getColStride() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept;
    };

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ device_obj<TensorSlice<X, DimR, DimC>>::device_obj(Ref tensor_, IndexVar auto... indices) : tensor(asStruct(tensor_)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(indices < tensor_.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ void device_obj<TensorSlice<X, DimR, DimC>>::resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow() && col == getCol());
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ size_t device_obj<TensorSlice<X, DimR, DimC>>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return getRow();
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ constexpr size_t device_obj<TensorSlice<X, DimR, DimC>>::getRowStride() const noexcept {
        return tensor.getDerived().getStride(DimR);
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ constexpr size_t device_obj<TensorSlice<X, DimR, DimC>>::getColStride() const noexcept {
        return tensor.getDerived().getStride(DimC);
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ auto device_obj<TensorSlice<X, DimR, DimC>>::data_handle(this auto&& self) noexcept {
        auto idx = self.index;
        idx[DimR] = 0;
        idx[DimC] = 0;
        return self.tensor.getDerived().data_ptr(idx);
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ consteval size_t device_obj<TensorSlice<X, DimR, DimC>>::getRowStrideAtCompile() noexcept {
        return std::remove_cvref_t<X>::getStrideAtCompile()[DimR];
    }

    template<Tensor X, int DimR, int DimC> requires(std::remove_cvref_t<X>::isCompact())
    __host__ __device__ consteval size_t device_obj<TensorSlice<X, DimR, DimC>>::getColStrideAtCompile() noexcept {
        return std::remove_cvref_t<X>::getStrideAtCompile()[DimC];
    }
}
