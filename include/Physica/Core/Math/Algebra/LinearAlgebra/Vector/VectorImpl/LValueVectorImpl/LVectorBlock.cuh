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
    template<Vector V, size_t Length>
    class device_obj<LVectorBlock<V, Length>> : public device_obj<LValueVector<LVectorBlock<V, Length>>> {
        using host_obj = LVectorBlock<V, Length>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
        using DeviceVector = device_obj<V>;
    public:
        using ScalarType = Base::ScalarType;
    private:
        Physica::PlainStruct<DeviceVector> vec;
        size_t from;
        size_t to;
    public:
        __host__ __device__ device_obj(device_obj<LValueVector<V>>& vec_, size_t from_, size_t to_);
        __host__ __device__ device_obj(device_obj<LValueVector<V>>& vec_, size_t from_);
        device_obj(const This& block) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        This& operator=(const This& v);
        This& operator=(This&& v) noexcept;
        /* Operations */
        template<size_t Length_ = Dynamic>
        [[nodiscard]] __host__ __device__ auto head(this auto&&, size_t to = Length_) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] __host__ __device__ auto tail(this auto&&, size_t from) noexcept;
        template<size_t Length_ = Dynamic>
        [[nodiscard]] __host__ __device__ auto segment(this auto&&, size_t from, size_t to) noexcept;

        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t index) noexcept;
    };

    template<Vector V, size_t Length>
    __host__ __device__ device_obj<LVectorBlock<V, Length>>::device_obj(
            device_obj<LValueVector<V>>& vec_, size_t from_, size_t to_) : vec(asStruct(vec_)), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    __host__ __device__ device_obj<LVectorBlock<V, Length>>::device_obj(
            device_obj<LValueVector<V>>& vec_, size_t from_) : device_obj(vec_, from_, vec_.getLength()) {}

    template<Vector V, size_t Length>
    auto device_obj<LVectorBlock<V, Length>>::operator=(const This& v) -> This& {
        Base::operator=(static_cast<const RValueVector<This>&>(v));
        return *this;
    }

    template<Vector V, size_t Length>
    auto device_obj<LVectorBlock<V, Length>>::operator=(This&& v) noexcept -> This& {
        Base::operator=(static_cast<const RValueVector<This>&>(v));
        return *this;
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    __host__ __device__ auto device_obj<LVectorBlock<V, Length>>::head(this auto&& self, size_t to) noexcept {
        return device_obj<LVectorBlock<V, Length_>>(self.vec, self.from, self.from + to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    __host__ __device__ auto device_obj<LVectorBlock<V, Length>>::tail(this auto&& self, size_t from) noexcept {
        return device_obj<LVectorBlock<V, Length_>>(self.vec, self.from + from, self.to);
    }

    template<Vector V, size_t Length>
    template<size_t Length_>
    __host__ __device__ auto device_obj<LVectorBlock<V, Length>>::segment(this auto&& self, size_t from, size_t to) noexcept {
        return device_obj<LVectorBlock<V, Length_>>(self.vec, self.from + from, self.from + to);
    }

    template<Vector V, size_t Length>
    __host__ __device__ size_t device_obj<LVectorBlock<V, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }

    template<Vector V, size_t Length>
    __host__ __device__ auto device_obj<LVectorBlock<V, Length>>::data_ptr(this auto&& self, size_t index) noexcept {
        assert((self.from + index) < self.to);
        return self.vec.getDerived().data_ptr(self.from + index);
    }
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<device_obj<LVectorBlock<V, Length>>> : public Traits<LVectorBlock<V, Length>> {};
}
