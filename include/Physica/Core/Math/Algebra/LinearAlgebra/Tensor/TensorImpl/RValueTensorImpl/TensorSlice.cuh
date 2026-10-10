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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.cuh"
#include "../RValueTensor.cuh"
#include "TensorSlice.h"

namespace Physica {
    template<class X, int DimR, int DimC>
    class device_obj<TensorSlice<X, DimR, DimC>> : public device_obj<RValueMatrix<TensorSlice<X, DimR, DimC>>> {
        static_assert(DimR != DimC, "[Error]: DimR and DimC must be different");

        using host_obj = TensorSlice<X, DimR, DimC>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj_t<X>;
        using IndexType = std::remove_cvref_t<X>::IndexType;
    protected:
        using typename Base::T;
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
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return tensor.getDerived().dim(DimR); }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return tensor.getDerived().dim(DimC); }
        [[nodiscard]] __host__ __device__ size_t getOrder() const noexcept;
    };

    template<class X, int DimR, int DimC>
    __host__ __device__ device_obj<TensorSlice<X, DimR, DimC>>::device_obj(Ref tensor_, IndexVar auto... indices) : tensor(asStruct(tensor_)) {
        size_t i = 0;
        ([&]() {
            if constexpr (std::integral<decltype(indices)>) {
                assert(static_cast<size_t>(indices) < tensor_.dim(i));
                index[i] = indices;
            }
            i += 1;
        }(), ...);
    }

    template<class X, int DimR, int DimC>
    __device__ auto device_obj<TensorSlice<X, DimR, DimC>>::calc(size_t row, size_t col, [[maybe_unused]] instanceof_x<ThreadBlock> auto block) const -> T {
        assert(row < getRow());
        assert(col < getCol());
        auto idx = index;
        idx[DimR] = row;
        idx[DimC] = col;
        return tensor.getDerived().calc(idx);
    }

    template<class X, int DimR, int DimC>
    __host__ __device__ size_t device_obj<TensorSlice<X, DimR, DimC>>::getOrder() const noexcept {
        assert(Base::isSquare() && "[Error]: getOrder() assumes square matrix");
        return getRow();
    }
}

namespace Physica {
    template<class X, int DimR, int DimC>
    class Traits<device_obj<TensorSlice<X, DimR, DimC>>> : public Traits<TensorSlice<X, DimR, DimC>> {};
}
