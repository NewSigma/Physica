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

#include "../RValueMatrix.cuh"

namespace Physica {
    template<Matrix M>
    class device_obj<MinorDiagR<M>> final : public device_obj<RValueVector<MinorDiagR<M>>> {
        using host_obj = MinorDiagR<M>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        ssize_t shift;
    public:
        __host__ __device__ device_obj(Ref mat, ssize_t shift);
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operations */
        [[nodiscard]] __device__ decltype(auto) calc(size_t index) const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ auto&& getExpr(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getExpr().getRow() - std::abs(shift); }
    };

    template<Matrix M>
    __host__ __device__ device_obj<MinorDiagR<M>>::device_obj(Ref mat, ssize_t shift) : mat(asStruct(mat)), shift(shift) {
        assert(getExpr().isSquare());
        assert(std::abs(shift) < getExpr().getCol());
    }

    template<Matrix M>
    __device__ decltype(auto) device_obj<MinorDiagR<M>>::calc(size_t index) const noexcept {
        size_t r = shift < 0 ? -shift : 0;
        size_t c = shift > 0 ? shift : 0;
        return getExpr().calc(r + index, c + index);
    }

    template<Matrix M>
    __host__ __device__ auto&& device_obj<MinorDiagR<M>>::getExpr(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), M>(self.mat.getDerived());
    }
}

namespace Physica {
    template<Matrix M>
    class Traits<device_obj<MinorDiagR<M>>> : public Traits<MinorDiagR<M>> {};
}
