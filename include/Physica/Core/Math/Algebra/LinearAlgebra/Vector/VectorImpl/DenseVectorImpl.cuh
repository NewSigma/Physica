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
    __host__ __device__ device_obj<DenseVector<T, Length, Allocator>>::device_obj(size_t length) : Storage(length) {}

    template<Scalar T, size_t Length, class Allocator>
    __host__ __device__ device_obj<DenseVector<T, Length, Allocator>>::device_obj(size_t length, T init) : device_obj(length) {
        *this = init;
    }

    template<Scalar T, size_t Length, class Allocator>
    device_obj<DenseVector<T, Length, Allocator>>::device_obj(const host_obj& obj) : Storage(obj) {}

    template<Scalar T, size_t Length, class Allocator>
    __host__ __device__ device_obj<DenseVector<T, Length, Allocator>>::device_obj(const Vector auto& v) : Storage(v.getLength()) {
        v.assign(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::resize(const Vector auto& x) {
        resize(x.getLength());
    }

    template<Scalar T, size_t Length, class Allocator>
    auto device_obj<DenseVector<T, Length, Allocator>>::toHost() const -> host_obj {
        return host_obj(Storage::toHost());
    }

    template<Scalar T, size_t Length, class Allocator>
    auto device_obj<DenseVector<T, Length, Allocator>>::toHostAsync() const -> host_obj {
        return host_obj(Storage::toHostAsync());
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::toHost(host_obj& obj) const {
        Storage::toHost(obj);
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::toHostAsync(host_obj& obj) const {
        Storage::toHostAsync(obj);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    void device_obj<DenseVector<T, Length, Allocator>>::random_uniform() {
        host_obj::template random_uniform<R>(getLength()).toDevice(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    void device_obj<DenseVector<T, Length, Allocator>>::random_normal() {
        host_obj::template random_normal<R>(getLength()).toDevice(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    void device_obj<DenseVector<T, Length, Allocator>>::random_any(auto& distribution) {
        host_obj::template random_any<R>(getLength(), distribution).toDevice(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Storage::swap(obj);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    auto device_obj<DenseVector<T, Length, Allocator>>::random_uniform(size_t len) -> This {
        return host_obj::template random_uniform<R>(len).toDevice();
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    auto device_obj<DenseVector<T, Length, Allocator>>::random_normal(size_t len) -> This {
        return host_obj::template random_normal<R>(len).toDevice();
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    auto device_obj<DenseVector<T, Length, Allocator>>::random_any(size_t len, auto& distribution) -> This {
        return host_obj::template random_any<R>(len, distribution).toDevice();
    }
}
