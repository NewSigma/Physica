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

#include "../RValueMatrix.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<DiagVector<M>> : public device_obj<RValueVector<DiagVector<M>>> {
        using host_obj = DiagVector<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
    public:
        __host__ __device__ explicit device_obj(Ref mat) : mat(asStruct(mat)) {}
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Operations */
        using Base::calc;
        [[nodiscard]] __device__ T calc(size_t index, instanceof_x<ThreadBlock> auto block) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
    };

    template<Matrix M>
    __device__ auto device_obj<DiagVector<M>>::calc(size_t index, instanceof_x<ThreadBlock> auto block) const -> T {
        return getExpr().calc(index, index, block);
    }

    template<Matrix M>
    __host__ __device__ size_t device_obj<DiagVector<M>>::getLength() const noexcept {
        return std::min(getExpr().getCol(), getExpr().getRow());
    }

    template<Matrix M>
    __host__ __device__ auto&& device_obj<DiagVector<M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<DiagVector<M>>> : public Traits<DiagVector<M>> {
        static_assert(std::is_reference<M>::value);
    };
}
