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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/CompactVector.cuh"

namespace Physica {
    template<class V> requires(std::remove_cvref_t<V>::isCompact())
    class device_obj<RealVector<V>> : public device_obj<StridedVector<RealVector<V>>> {
        using host_obj = RealVector<V>;
        using This = device_obj<host_obj>;
        using Base = device_obj<StridedVector<host_obj>>;
        using Ref = add_device_obj_t<V>;
    protected:
        using typename Base::T;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<V>>> v;
    public:
        __host__ __device__ explicit device_obj(Ref v_) : v(asStruct(v_)) {}
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data_handle(this auto&& self) noexcept;
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return v.getDerived().getLength(); }
    };

    template<class V> requires(std::remove_cvref_t<V>::isCompact())
    __host__ __device__ auto device_obj<RealVector<V>>::data_handle(this auto&& self) noexcept {
        return self.v.getDerived()[0].real_ptr();
    }
}
