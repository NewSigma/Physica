/*
 * Copyright 2024-2025 Weibo He.
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

#include "RVectorBlock.h"

namespace Physica {
    template<Vector V, size_t Length>
    class device_obj<RVectorBlock<V, Length>> : public device_obj<RValueVector<RVectorBlock<V, Length>>> {
        using host_obj = RVectorBlock<V, Length>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using ScalarType = Base::ScalarType;
    private:
        const V& vec;
        size_t from;
        size_t to;
    public:
        __host__ __device__ device_obj(const device_obj<V>& vec_, size_t from_, size_t to_);
        __host__ __device__ device_obj(const device_obj<V>& vec_, size_t from_);
        device_obj(const This& block) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { return vec.calc(index + from); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
    };

    template<Vector V, size_t Length>
    __host__ __device__ device_obj<RVectorBlock<V, Length>>::device_obj(
            const device_obj<V>& vec_, size_t from_, size_t to_) : vec(asStruct(vec_)), from(from_), to(to_) {
        assert(from_ < to);
        assert(to <= vec.getLength());
        assert(Length == Dynamic || Length == getLength());
    }

    template<Vector V, size_t Length>
    __host__ __device__ device_obj<RVectorBlock<V, Length>>::device_obj(
            const device_obj<V>& vec_, size_t from_) : device_obj(vec_, from_, vec_.getLength()) {}

    template<Vector V, size_t Length>
    __host__ __device__ size_t device_obj<RVectorBlock<V, Length>>::getLength() const noexcept {
        if constexpr (Length == Dynamic)
            return to - from;
        return Length;
    }
}

namespace Physica {
    template<Vector V, size_t Length>
    class Traits<device_obj<RVectorBlock<V, Length>>> : public Traits<RVectorBlock<V, Length>> {};
}
