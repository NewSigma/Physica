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

#include "../LValueMatrix.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<MinorDiagL<M>> final : public device_obj<LValueVector<MinorDiagL<M>>> {
        using host_obj = MinorDiagL<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        ssize_t shift;
    public:
        __host__ __device__ device_obj(Ref mat, ssize_t shift);
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t size) { assert(getLength() == size); }
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t index) noexcept;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getExpr().getRow() - std::abs(shift); }
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&&) noexcept;
    };

    template<Matrix M>
    __host__ __device__ device_obj<MinorDiagL<M>>::device_obj(Ref mat, ssize_t shift) : mat(asStruct(mat)), shift(shift) {
        assert(getExpr().isSquare());
        assert(std::abs(shift) < getExpr().getRow());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<MinorDiagL<M>>::data_ptr(this auto&& self, size_t index) noexcept {
        auto shift = self.shift;
        size_t r = shift < 0 ? -shift : 0;
        size_t c = shift > 0 ? shift : 0;
        return self.getExpr().data_ptr(r + index, c + index);
    }

    template<Matrix M>
    __host__ __device__ auto&& device_obj<MinorDiagL<M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref>(self.mat.getDerived());
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<MinorDiagL<M>>> : public Traits<MinorDiagL<M>> {};
}
