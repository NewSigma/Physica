/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/PlainStruct.h"
#include "GEVM.h"

namespace Physica {
    template<Vector V, Matrix M>
    class device_obj<GEVM<V, M>>
            : public device_obj<RValueMatrix<GEVM<V, M>>> {
        static_assert(M::RowAtCompile == 1 || M::RowAtCompile == Dynamic,
                      "Row and column do not match in matrix product");
        using host_obj = GEVM<V, M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    protected:
        using typename Base::T;
    private:
        Physica::PlainStruct<const device_obj<V>> vec;
        Physica::PlainStruct<const device_obj<M>> mat;
    public:
        __host__ __device__ device_obj(const device_obj<V>& vec_, const device_obj<M>& mat_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const { return vec.getDerived().getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getDerived().getCol(); }
    };

    template<Vector V, Matrix M>
    __host__ __device__ device_obj<GEVM<V, M>>::device_obj(const device_obj<V>& vec_, const device_obj<M>& mat_)
            : vec(asStruct(vec_)), mat(asStruct(mat_)) {
        assert(mat_.getRow() == 1);
    }

    template<Vector V, Matrix M>
    __device__ auto device_obj<GEVM<V, M>>::calc(size_t row, size_t col) const -> T {
        return vec.getDerived().calc(row) * mat.getDerived().calc(0, col);
    }

    template<Vector V, Matrix M>
    [[nodiscard]] __host__ __device__ auto operator*(const device_obj<V>& vec, const device_obj<M>& mat) noexcept requires(M::RowAtCompile == 1) {
        assert(mat.getRow() == 1);
        return device_obj<GEVM<V, M>>(vec, mat);
    }
}

namespace Physica {
    template<Vector V, Matrix M>
    class Traits<device_obj<GEVM<V, M>>> : public Traits<GEVM<V, M>> {};
}
