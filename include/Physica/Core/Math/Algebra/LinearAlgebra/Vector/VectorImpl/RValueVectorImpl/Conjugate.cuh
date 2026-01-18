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

#include "../RValueVector.cuh"

namespace Physica {
    template<Vector V>
    class device_obj<Conjugate<V>> : public device_obj<RValueVector<Conjugate<V>>> {
        using host_obj = Conjugate<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj<V>::type;
    protected:
        using typename Base::T;
        using typename Base::Tv;
        using typename Base::Tm;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> vec;
    public:
        __host__ __device__ explicit device_obj(Ref vec_) : vec(asStruct(vec_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t index) const { return getExpr().calc(index).conjugate(); }
        [[nodiscard]] __device__ Tv calc_value(size_t index) const { return getExpr().calc_value(index).conjugate(); }

        [[nodiscard]] __host__ __device__ decltype(auto) conjugate(this auto&&) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ const auto& getExpr() const noexcept { return vec.getDerived(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return getExpr().getLength(); }
    };

    template<Vector V>
    __host__ __device__ decltype(auto) device_obj<Conjugate<V>>::conjugate(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.vec);
    }
}

namespace Physica {
    template<Vector V>
    class Traits<device_obj<Conjugate<V>>> : public Traits<Conjugate<V>> {
        static_assert(!is_device_obj<V>::value);
    };
}

#ifdef PHYSICA_MKL
    #include "Conjugate_MKL.h"
#endif
