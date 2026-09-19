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

#include "../DenseTensor.cuh"

namespace Physica {
    template<Scalar T, int... Dims>
    device_obj<DenseTensor<T, Dims...>>::device_obj(Storage storage) noexcept : storage(std::move(storage)) {}

    template<Scalar T, int... Dims>
    __host__ __device__ device_obj<DenseTensor<T, Dims...>>::device_obj(IndexType shape, auto&&... args) {
        resize(std::move(shape), std::forward<decltype(args)>(args)...);
    }

    template<Scalar T, int... Dims>
    __host__ __device__ device_obj<DenseTensor<T, Dims...>>::device_obj(std::integral auto... dims) : storage(dims...) {}

    template<Scalar T, int... Dims>
    device_obj<DenseTensor<T, Dims...>>::device_obj(const Tensor auto& x) : This(x.getShape()) {
        x.assign(*this);
    }

    template<Scalar T, int... Dims>
    device_obj<DenseTensor<T, Dims...>>::device_obj(const host_obj& obj) : storage(obj.storage) {}

    template<Scalar T, int... Dims>
    __host__ __device__ void device_obj<DenseTensor<T, Dims...>>::resize(this auto& self, IndexType shape) {
        self.storage.resize(std::move(shape));
    }

    template<Scalar T, int... Dims>
    void device_obj<DenseTensor<T, Dims...>>::reserve(this auto& self, size_t size) {
        self.storage.reserve(size);
    }

    template<Scalar T, int... Dims>
    auto device_obj<DenseTensor<T, Dims...>>::toHost() const -> host_obj {
        return host_obj(storage.toHost());
    }

    template<Scalar T, int... Dims>
    auto device_obj<DenseTensor<T, Dims...>>::toHostAsync() const -> host_obj {
        return host_obj(storage.toHostAsync());
    }

    template<Scalar T, int... Dims>
    void device_obj<DenseTensor<T, Dims...>>::toHost(host_obj& obj) const {
        storage.toHost(obj.storage);
    }

    template<Scalar T, int... Dims>
    void device_obj<DenseTensor<T, Dims...>>::toHostAsync(host_obj& obj) const {
        storage.toHostAsync(obj.storage);
    }

    template<Scalar T, int... Dims>
    void device_obj<DenseTensor<T, Dims...>>::zeros(this auto& self) noexcept {
        self.storage.zeros();
    }

    template<Scalar T, int... Dims>
    void device_obj<DenseTensor<T, Dims...>>::junk(this auto& self) noexcept {
        self.storage.junk();
    }

    template<Scalar T, int... Dims>
    void device_obj<DenseTensor<T, Dims...>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        storage.swap(obj.storage);
    }

    template<Scalar T, int... Dims>
    template<RNG R>
    auto device_obj<DenseTensor<T, Dims...>>::random_uniform(IndexType shape) -> This {
        return host_obj::template random_uniform<R>(std::move(shape)).toDevice();
    }

    template<Scalar T, int... Dims>
    template<RNG R>
    auto device_obj<DenseTensor<T, Dims...>>::random_normal(IndexType shape) -> This {
        return host_obj::template random_normal<R>(std::move(shape)).toDevice();
    }
}
