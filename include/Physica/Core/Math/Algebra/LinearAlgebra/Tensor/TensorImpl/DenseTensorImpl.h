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
    DenseTensor<T, Dims...>::DenseTensor(std::integral auto... dims) : storage(dims...) {}

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
}
