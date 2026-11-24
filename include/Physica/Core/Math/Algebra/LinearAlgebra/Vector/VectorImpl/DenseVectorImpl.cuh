/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../DenseVector.cuh"

namespace Physica {
    template<Scalar T, size_t Length, class Allocator>
    device_obj<DenseVector<T, Length, Allocator>>::device_obj(Storage storage) noexcept : storage(std::move(storage)) {}

    template<Scalar T, size_t Length, class Allocator>
    __host__ __device__ device_obj<DenseVector<T, Length, Allocator>>::device_obj(size_t length) : storage(length) {}

    template<Scalar T, size_t Length, class Allocator>
    __host__ __device__ device_obj<DenseVector<T, Length, Allocator>>::device_obj(size_t length, T init) : device_obj(length) {
        *this = init;
    }

    template<Scalar T, size_t Length, class Allocator>
    device_obj<DenseVector<T, Length, Allocator>>::device_obj(const host_obj& obj) : storage(obj.storage) {}

    template<Scalar T, size_t Length, class Allocator>
    __host__ __device__ device_obj<DenseVector<T, Length, Allocator>>::device_obj(const Vector auto& v) : storage(v.getLength()) {
        v.assign(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    __host__ __device__ void device_obj<DenseVector<T, Length, Allocator>>::resize(this auto& self, const Vector auto& x) {
        self.resize(x.getLength());
    }

    template<Scalar T, size_t Length, class Allocator>
    __host__ __device__ void device_obj<DenseVector<T, Length, Allocator>>::resize(this auto& self, size_t size, auto&&... args) noexcept {
        self.storage.resize(size, std::forward<decltype(args)>(args)...);
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::reserve(this auto& self, size_t size) noexcept {
        self.storage.reserve(size);
    }

    template<Scalar T, size_t Length, class Allocator>
    auto device_obj<DenseVector<T, Length, Allocator>>::toHost() const -> host_obj {
        return host_obj(storage.toHost());
    }

    template<Scalar T, size_t Length, class Allocator>
    auto device_obj<DenseVector<T, Length, Allocator>>::toHostAsync() const -> host_obj {
        return host_obj(storage.toHostAsync());
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::toHost(host_obj& obj) const {
        storage.toHost(obj.storage);
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::toHostAsync(host_obj& obj) const {
        storage.toHostAsync(obj.storage);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    void device_obj<DenseVector<T, Length, Allocator>>::random_uniform(this auto& self) {
        host_obj::template random_uniform<R>(self.getLength()).toDevice(self);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    void device_obj<DenseVector<T, Length, Allocator>>::random_normal(this auto& self) {
        host_obj::template random_normal<R>(self.getLength()).toDevice(self);
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    void device_obj<DenseVector<T, Length, Allocator>>::random_any(this auto& self, auto& distribution) {
        host_obj::template random_any<R>(self.getLength(), distribution).toDevice(self);
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::zeros(this auto& self) noexcept {
        self.storage.zeros();
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::junk(this auto& self) noexcept {
        self.storage.junk();
    }

    template<Scalar T, size_t Length, class Allocator>
    void device_obj<DenseVector<T, Length, Allocator>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        storage.swap(obj.storage);
    }

    template<Scalar T, size_t Length, class Allocator>
    __host__ __device__ auto* device_obj<DenseVector<T, Length, Allocator>>::data_handle(this auto&& self) noexcept {
        return self.storage.data();
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
