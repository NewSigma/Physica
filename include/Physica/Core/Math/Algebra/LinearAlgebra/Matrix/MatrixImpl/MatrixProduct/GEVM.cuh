/*
 * Copyright 2023-2026 Weibo He.
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

#include "Physica/Core/Parallel/ThreadBlock.cuh"
#include "Physica/PlainStruct.h"
#include "GEVM.h"

namespace Physica {
    template<Vector V, Matrix M>
    class device_obj<GEVM<V, M>> : public device_obj<RValueMatrix<GEVM<V, M>>> {
        using host_obj = GEVM<V, M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    protected:
        using typename Base::T;
    private:
        Physica::PlainStruct<const std::remove_cvref_t<V>> vec;
        Physica::PlainStruct<const std::remove_cvref_t<M>> mat;
    public:
        __host__ __device__ device_obj(V vec_, M mat_);
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
        [[nodiscard]] __host__ __device__ size_t getRow() const { return vec.getDerived().getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat.getDerived().getCol(); }
    };

    template<Vector V, Matrix M>
    __host__ __device__ device_obj<GEVM<V, M>>::device_obj(V vec_, M mat_) : vec(asStruct(vec_)), mat(asStruct(mat_)) {}

    template<Vector V, Matrix M>
    __device__ auto device_obj<GEVM<V, M>>::calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto block) const -> T {
        return vec.getDerived().calc(row, block) * mat.getDerived().calc(0, col, block);
    }
}

namespace Physica {
    template<Vector V, Matrix M>
    class Traits<device_obj<GEVM<V, M>>> : public Traits<GEVM<V, M>> {};
}
