/*
 * Copyright 2024-2026 Weibo He.
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

#include "../LValueVector.cuh"

namespace Physica {
    template<Vector T, size_t Length>
    class device_obj<LVectorBlock<T, Length>> : public device_obj<LValueVector<LVectorBlock<T, Length>>> {
        using host_obj = LVectorBlock<T, Length>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
        using DeviceVector = device_obj<T>;
    public:
        using ScalarType = Base::ScalarType;
    private:
        Physica::PlainStruct<DeviceVector> vec;
        size_t from;
        size_t to;
    public:
        __host__ __device__ device_obj(device_obj<LValueVector<T>>& vec_, size_t from_, size_t to_);
        __host__ __device__ device_obj(device_obj<LValueVector<T>>& vec_, size_t from_);
        device_obj(const This& block) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& v);
        This& operator=(This&& v) noexcept;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t index) noexcept;
    };

    template<Vector T, size_t Length>
    __host__ __device__ device_obj<LVectorBlock<T, Length>>::device_obj(
            device_obj<LValueVector<T>>& vec_, size_t from_, size_t to_) : vec(asStruct(vec_)), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector T, size_t Length>
    __host__ __device__ device_obj<LVectorBlock<T, Length>>::device_obj(
            device_obj<LValueVector<T>>& vec_, size_t from_) : device_obj(vec_, from_, vec_.getLength()) {}

    template<Vector T, size_t Length>
    auto device_obj<LVectorBlock<T, Length>>::operator=(const This& v) -> This& {
        Base::operator=(static_cast<const RValueVector<This>&>(v));
        return *this;
    }

    template<Vector T, size_t Length>
    auto device_obj<LVectorBlock<T, Length>>::operator=(This&& v) noexcept -> This& {
        Base::operator=(static_cast<const RValueVector<This>&>(v));
        return *this;
    }

    template<Vector T, size_t Length>
    __host__ __device__ size_t device_obj<LVectorBlock<T, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }

    template<Vector T, size_t Length>
    __host__ __device__ auto device_obj<LVectorBlock<T, Length>>::data_ptr(this auto&& self, size_t index) noexcept {
        assert((self.from + index) < self.to);
        return self.vec.getDerived().data_ptr(self.from + index);
    }
}

namespace Physica {
    template<Vector T, size_t Length>
    class Traits<device_obj<LVectorBlock<T, Length>>> : public Traits<LVectorBlock<T, Length>> {};
}
