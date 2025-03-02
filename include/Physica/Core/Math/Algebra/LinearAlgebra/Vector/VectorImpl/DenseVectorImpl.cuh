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

#include "../DenseVector.cuh"

namespace Physica {
    template<Scalar T, size_t Length, class Allocator>
    device_obj<DenseVector<T, Length, Allocator>>::device_obj(const host_obj& obj) : Storage(obj) {}

    template<Scalar T, size_t Length, class Allocator>
    template<Vector V>
    __host__ __device__ device_obj<DenseVector<T, Length, Allocator>>::device_obj(const V& v) requires(CUDA<V>) : Storage(v.getLength()) {
        v.assign(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<Vector V>
    void device_obj<DenseVector<T, Length, Allocator>>::resize(const V& x) {
        resize(x.getLength());
    }

    template<Scalar T, size_t Length, class Allocator>
    inline auto device_obj<DenseVector<T, Length, Allocator>>::toHost() const -> host_obj {
        return host_obj(Storage::toHost());
    }

    template<Scalar T, size_t Length, class Allocator>
    inline auto device_obj<DenseVector<T, Length, Allocator>>::toHostAsync() const -> host_obj {
        return host_obj(Storage::toHostAsync());
    }

    template<Scalar T, size_t Length, class Allocator>
    inline void device_obj<DenseVector<T, Length, Allocator>>::toHost(host_obj& obj) const {
        Storage::toHost(obj);
    }

    template<Scalar T, size_t Length, class Allocator>
    inline void device_obj<DenseVector<T, Length, Allocator>>::toHostAsync(host_obj& obj) const {
        Storage::toHostAsync(obj);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    inline void device_obj<DenseVector<T, Length, Allocator>>::random_uniform() {
        host_obj::template random_uniform<R>(getLength()).toDevice(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    inline void device_obj<DenseVector<T, Length, Allocator>>::random_normal() {
        host_obj::template random_normal<R>(getLength()).toDevice(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R, class Distribution>
    inline void device_obj<DenseVector<T, Length, Allocator>>::random_any(Distribution& dist) {
        host_obj::template random_any<R, Distribution>(getLength(), dist).toDevice(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Storage::swap(obj);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    inline auto device_obj<DenseVector<T, Length, Allocator>>::random_uniform(size_t len) -> This {
        return host_obj::template random_uniform<R>(len).toDevice();
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    inline auto device_obj<DenseVector<T, Length, Allocator>>::random_normal(size_t len) -> This {
        return host_obj::template random_normal<R>(len).toDevice();
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R, class Distribution>
    inline auto device_obj<DenseVector<T, Length, Allocator>>::random_any(size_t len, Distribution& dist) -> This {
        return host_obj::template random_any<R, Distribution>(len, dist).toDevice();
    }
}
