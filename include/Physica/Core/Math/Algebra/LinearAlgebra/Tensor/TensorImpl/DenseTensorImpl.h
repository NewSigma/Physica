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

#include "../DenseTensor.h"

namespace Physica {
    template<Scalar T, int... Dims>
    DenseTensor<T, Dims...>::DenseTensor(ArrayND<T, Dims...> storage) noexcept : storage(std::move(storage)) {}

    template<Scalar T, int... Dims>
    DenseTensor<T, Dims...>::DenseTensor(IndexType shape, auto&&... args) : storage(std::move(shape), std::forward<decltype(args)>(args)...) {}

    template<Scalar T, int... Dims>
    DenseTensor<T, Dims...>::DenseTensor(std::same_as<size_t> auto... dims) : storage(dims...) {}

    template<Scalar T, int... Dims>
    DenseTensor<T, Dims...>::DenseTensor(const Tensor auto& x) : This(x.getShape()) {
        x.assign(*this);
    }

    template<Scalar T, int... Dims>
    void DenseTensor<T, Dims...>::resize(this auto& self, IndexType shape) {
        self.storage.resize(std::move(shape));
    }

    template<Scalar T, int... Dims>
    void DenseTensor<T, Dims...>::reserve(this auto& self, size_t size) noexcept {
        self.storage.reserve(size);
    }

    template<Scalar T, int... Dims>
    void DenseTensor<T, Dims...>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        storage.swap(obj.storage);
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::data_handle(this auto&& self) noexcept {
        return self.storage.data();
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::getShape() const noexcept -> IndexType {
        return storage.getShape();
    }

    template<Scalar T, int... Dims>
    constexpr size_t DenseTensor<T, Dims...>::dim(int index) const noexcept {
        return storage.dim(index);
    }

    template<Scalar T, int... Dims>
    auto&& DenseTensor<T, Dims...>::asArray(this auto&& self) noexcept {
        return self.storage.asArray();
    }

    template<Scalar T, int... Dims>
    size_t DenseTensor<T, Dims...>::getSize() const noexcept {
        return storage.getSize();
    }

    template<Scalar T, int... Dims>
    __host__ __device__ consteval auto DenseTensor<T, Dims...>::getStrideAtCompile() noexcept -> IndexType {
        IndexType result{};
        if constexpr (sizeof...(Dims) == 1) {
            for (int i = 0; i < NDim - 1; ++i)
                result[i] = Dynamic;
            result[NDim - 1] = 1;
        }
        else {
            constexpr std::array<int, NDim> dims{Dims...};
            size_t stride = 1;
            bool known = true;
            for (int i = NDim - 1; i >= 0; --i) {
                result[i] = known ? stride : Dynamic;
                if (dims[i] == Dynamic)
                    known = false;
                else
                    stride *= static_cast<size_t>(dims[i]);
            }
        }
        return result;
    }

    template<Scalar T, int... Dims>
    __host__ __device__ consteval size_t DenseTensor<T, Dims...>::getSizeAtCompile() noexcept {
        return ArrayND<T, Dims...>::SizeAtCompile;
    }

    template<Scalar T, int... Dims>
    template<RNG R>
    auto DenseTensor<T, Dims...>::random_uniform(IndexType shape) -> This {
        auto result = This(std::move(shape));
        result.asArray().template random_uniform<R>();
        return result;
    }

    template<Scalar T, int... Dims>
    template<RNG R>
    auto DenseTensor<T, Dims...>::random_normal(IndexType shape) -> This {
        auto result = This(std::move(shape));
        result.asArray().template random_normal<R>();
        return result;
    }

    template<Scalar T, int... Dims>
    void DenseTensor<T, Dims...>::zeros(this auto& self) noexcept {
        self.storage.zeros();
    }

    template<Scalar T, int... Dims>
    void DenseTensor<T, Dims...>::junk(this auto& self) noexcept {
        self.storage.junk();
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::zeros(IndexType shape) -> This {
        auto result = This(std::move(shape));
        result.zeros();
        return result;
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::zeros(std::same_as<size_t> auto... dims) -> This {
        static_assert(sizeof...(dims) == ArrayND<T, Dims...>::NDim, "[Error]: NDim is not consistent");
        return zeros(IndexType({dims...}));
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::junk(IndexType shape) -> This {
        auto result = This(std::move(shape));
        result.junk();
        return result;
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::junk(std::same_as<size_t> auto... dims) -> This {
        static_assert(sizeof...(dims) == ArrayND<T, Dims...>::NDim, "[Error]: NDim is not consistent");
        return junk(IndexType({dims...}));
    }

    template<Scalar T, int... Dims>
    template<RNG R>
    auto DenseTensor<T, Dims...>::random_uniform(std::same_as<size_t> auto... dims) -> This {
        static_assert(sizeof...(dims) == ArrayND<T, Dims...>::NDim, "[Error]: NDim is not consistent");
        return random_uniform<R>(IndexType({dims...}));
    }

    template<Scalar T, int... Dims>
    template<RNG R>
    auto DenseTensor<T, Dims...>::random_normal(std::same_as<size_t> auto... dims) -> This {
        static_assert(sizeof...(dims) == ArrayND<T, Dims...>::NDim, "[Error]: NDim is not consistent");
        return random_normal<R>(IndexType({dims...}));
    }

    template<Scalar T, int... Dims>
    template<RNG R>
    auto DenseTensor<T, Dims...>::random_any(IndexType shape, auto& distribution) -> This {
        auto result = This(std::move(shape));
        result.asArray().template random_any<R>(distribution);
        return result;
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::generate(std::invocable<IndexType> auto fn, IndexType shape) -> This {
        auto result = This(std::move(shape));
        result.forND([&fn](T& value, IndexType index) { value = fn(index); });
        return result;
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::generate(std::invocable<IndexType> auto fn, std::same_as<size_t> auto... dims) -> This {
        static_assert(sizeof...(dims) == ArrayND<T, Dims...>::NDim, "[Error]: NDim is not consistent");
        return generate(std::move(fn), IndexType({dims...}));
    }

    template<Scalar T, int... Dims>
    auto DenseTensor<T, Dims...>::read(IndexType shape, const T* __restrict p) noexcept -> This {
        return This(ArrayND<T, Dims...>::read(std::move(shape), p));
    }
}
