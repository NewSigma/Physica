/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.cuh"

namespace Physica {
    template<class V>
    class device_obj<RealVectorR<V>> : public device_obj<RValueVector<RealVectorR<V>>> {
        using host_obj = RealVectorR<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<V>> v;
    public:
        __host__ __device__ explicit device_obj(const device_obj<V>& v_) : v(asStruct(v_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t s) const { return v.getDerived().calc(s).real(); }
        [[nodiscard]] __device__ Tv calc_value(size_t s) const { return v.getDerived().calc_value(s).real(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };
}

namespace Physica {
    template<class V>
    class Traits<device_obj<RealVectorR<V>>> : public Traits<RealVectorR<V>> {};
}
