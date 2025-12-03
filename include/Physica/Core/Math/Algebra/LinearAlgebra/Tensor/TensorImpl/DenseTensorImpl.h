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
    DenseTensor<T, Dims...>::DenseTensor(std::integral auto... dims) : Storage(dims...) {
        static_assert(sizeof...(dims) == NDim, "[Error]: NDim is not consistent");
    }

    template<Scalar T, int... Dims>
    DenseTensor<T, Dims...>::DenseTensor(IndexType shape, auto&&... args) : Storage(std::move(shape), std::forward<decltype(args)>(args)...) {}

    template<Scalar T, int... Dims>
    void DenseTensor<T, Dims...>::resize(IndexType shape) {
        Storage::resize(std::move(shape));
    }

    template<Scalar T, int... Dims>
    void DenseTensor<T, Dims...>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Storage::swap(obj);
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
