/*
 * Copyright 2023-2025 Weibo He.
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

#include "../LValueTensor.h"

namespace Physica {
    template<class Derived>
    template<Scalar U>
    Derived& LValueTensor<Derived>::operator=(const U& x) requires(!isReverseDiff || !ReverseDiff<U>) {
        flatten() = x;
        return Base::getDerived();
    }

    template<class Derived>
    Derived& LValueTensor<Derived>::operator=(const Tensor auto& other) {
        if constexpr (std::is_same<Derived, std::remove_cvref_t<decltype(other)>>::value)
            assert(this != &other && "[Error]: Self assign is likely a bug");
        Derived& x = Base::getDerived();
        x.resize(other);
        other.assign(x);
        return x;
    }

    template<class Derived>
    void LValueTensor<Derived>::operator+=(const Tensor auto& x) {
        Base::getDerived() = Base::getDerived() + x;
    }

    template<class Derived>
    void LValueTensor<Derived>::operator-=(const Tensor auto& x) {
        Base::getDerived() = Base::getDerived() - x;
    }

    template<class Derived>
    void LValueTensor<Derived>::forND(std::invocable<T&, IndexType> auto fn) {
        Physica::forND(Base::getShape(), [this, fn](const IndexType& index) {
            fn(operator()(index), index);
        });
    }

    template<class Derived>
    void LValueTensor<Derived>::forND(std::invocable<const T&, IndexType> auto fn) const {
        Physica::forND(Base::getShape(), [this, fn](const IndexType& index) {
            fn(operator()(index), index);
        });
    }

    template<class Derived>
    LTensorBlock<Derived> LValueTensor<Derived>::block(Index3D from, Index3D count) {
        return {Base::getDerived(), from, count};
    }

    template<class Derived>
    const LTensorBlock<Derived> LValueTensor<Derived>::block(Index3D from, Index3D count) const {
        return {Base::getDerived(), from, count};
    }

    template<class Derived>
    auto LValueTensor<Derived>::flatten() {
        return FlattenL<Derived>(Base::getDerived());
    }

    template<class Derived>
    const auto LValueTensor<Derived>::flatten() const {
        return FlattenL<Derived>(Base::getConstCastDerived());
    }

    template<class Derived>
    template<RNG R>
    void LValueTensor<Derived>::random_uniform() {
        flatten().template random_uniform<R>();
    }

    template<class Derived>
    template<RNG R>
    void LValueTensor<Derived>::random_normal() {
        flatten().template random_normal<R>();
    }

    template<class Derived>
    void operator*=(LValueTensor<Derived>& grid, typename Derived::ScalarType factor) {
        for (size_t i = 0; i < grid.getDimX(); ++i)
            for (size_t j = 0; j < grid.getDimY(); ++j)
                for (size_t k = 0; k < grid.getDimZ(); ++k)
                    grid(i, j, k) *= factor;
    }
}
