/*
 * Copyright 2023-2026 Weibo He.
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
#include "TensorFiber.h"
#include "TensorSlice.h"

namespace Physica {
    template<class Derived>
    template<Scalar U>
    Derived& LValueTensor<Derived>::operator=(const U& x) requires(!isReverseDiff || !ReverseDiff<U>) {
        flatten() = x;
        return Base::getDerived();
    }

    template<class Derived>
    void LValueTensor<Derived>::operator+=(const Scalar auto& x) {
        auto& t = Base::getDerived();
        (t + x).assign(t);
    }

    template<class Derived>
    void LValueTensor<Derived>::operator-=(const Scalar auto& x) {
        auto& t = Base::getDerived();
        (t - x).assign(t);
    }

    template<class Derived>
    void LValueTensor<Derived>::operator*=(const Scalar auto& x) {
        auto& t = Base::getDerived();
        (t * x).assign(t);
    }

    template<class Derived>
    void LValueTensor<Derived>::operator/=(const Scalar auto& x) {
        auto& t = Base::getDerived();
        (t / x).assign(t);
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
    decltype(auto) LValueTensor<Derived>::operator[](this auto&& self, size_t x, size_t y, size_t z) {
        return *self.data_ptr({x, y, z});
    }

    template<class Derived>
    decltype(auto) LValueTensor<Derived>::operator[](this auto&& self, Index3D index) {
        return *self.data_ptr(index);
    }

    template<class Derived>
    void LValueTensor<Derived>::forND(std::invocable<T&, IndexType> auto fn) {
        Physica::forND(Base::getShape(), [this, fn](const IndexType& index) {
            fn(operator[](index), index);
        });
    }

    template<class Derived>
    void LValueTensor<Derived>::forND(std::invocable<const T&, IndexType> auto fn) const {
        Physica::forND(Base::getShape(), [this, fn](const IndexType& index) {
            fn(operator[](index), index);
        });
    }

    template<class Derived>
    auto LValueTensor<Derived>::fiber(this auto&& self, int dim, IndexType index) noexcept {
        using Self = decltype(self);
        return TensorFiber<Self>(std::forward<Self>(self), dim, index);
    }

    template<class Derived>
    auto LValueTensor<Derived>::slice(this auto&& self, int dimRow, int dimCol, IndexType index) noexcept {
        using Self = decltype(self);
        return TensorSlice<Self>(std::forward<Self>(self), dimRow, dimCol, index);
    }

    template<class Derived>
    auto LValueTensor<Derived>::block(this auto&& self, Index3D from, Index3D count) noexcept {
        using Self = decltype(self);
        return LTensorBlock<Self>(std::forward<Self>(self), from, count);
    }

    template<class Derived>
    auto LValueTensor<Derived>::flatten(this auto&& self) {
        using Self = decltype(self);
        return FlattenL<Self>(std::forward<Self>(self));
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
    auto LValueTensor<Derived>::data_ptr(this auto&& self, Index3D index) noexcept {
        return self.getDerived().data_ptr(index);
    }

    template<class Derived>
    void operator*=(LValueTensor<Derived>& grid, typename Derived::ScalarType factor) {
        for (size_t i = 0; i < grid.getDimX(); ++i)
            for (size_t j = 0; j < grid.getDimY(); ++j)
                for (size_t k = 0; k < grid.getDimZ(); ++k)
                    grid(i, j, k) *= factor;
    }
}
