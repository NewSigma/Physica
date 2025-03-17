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

#include "../RValueTensor.h"
#include "Physica/Core/Utils/Container/ArrayND.h"

namespace Physica {
    template<class Derived>
    template<Tensor X>
    void RValueTensor<Derived>::assign(X& x) const {
        const auto& derived = Base::getDerived();
        for (size_t i = 0; i < derived.getSize(); ++i) {
            const auto indices = toIndexND(i);
            x(indices) = calc(indices);
        }
    }

    template<class Derived>
    size_t RValueTensor<Derived>::toIndex1D(const IndexArray& indices) const noexcept {
        return ArrayND<ScalarType, Dim>::toIndex1D(getShape(), indices);
    }

    template<class Derived>
    auto RValueTensor<Derived>::toIndexND(size_t index) const noexcept -> IndexArray {
        return ArrayND<ScalarType, Dim>::toIndexND(getShape(), index);
    }

    template<class Derived>
    template<class Functor>
    void RValueTensor<Derived>::forND(Functor func) const {
        Physica::forND(getShape(), [this, func](const IndexArray& index) {
            func(calc(index), index);
        });
    }

    template<class Derived>
    auto RValueTensor<Derived>::reals() const noexcept {
        return RealTensor<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueTensor<Derived>::imags() const noexcept {
        return ImagTensor<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueTensor<Derived>::squaredNorms() const noexcept {
        return SquaredNormTensor<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueTensor<Derived>::norms() const noexcept {
        return NormTensor<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueTensor<Derived>::values() const noexcept {
        return ValueTensor<Derived>(Base::getDerived());
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueTensor<Derived>::grads() const noexcept {
        if constexpr (isReverseDiff)
            return Base::getDerived().template grads<GradOrder>();
        else
            return grads_impl<GradOrder>();
    }

    template<class Derived>
    auto RValueTensor<Derived>::getShape() const noexcept -> IndexArray {
        const auto& x = Base::getDerived();
        IndexArray shape(x.getDim());
        for (int i = 0; i < static_cast<int>(shape.getLength()); ++i)
            shape[i] = x.getShape(i);
        return shape;
    }

    template<class Derived>
    auto RValueTensor<Derived>::getDim() const noexcept {
        if constexpr (Dim == Dynamic)
            return Base::getDerived().getDim();
        else
            return Dim;
    }

    template<class Derived>
    size_t RValueTensor<Derived>::getSize() const noexcept {
        size_t size = getShape(0);
        for (int i = 1; i < getDim(); ++i)
            size *= getShape(i);
        return size;
    }

    template<class Derived>
    template<Scalar T, bool IsUnitLattice, class Functor>
    inline void RValueTensor<Derived>::forPointIndexInTensor(
            const RValueTensor& grid, const PeriodicCell<T, 3>::LatticeMatrix& lattice, Functor func) {
        forPointIndexInTensor<T, IsUnitLattice, Functor>(grid.getDim(), lattice, func);
    }

    template<class Derived>
    size_t RValueTensor<Derived>::toSize(const IndexArray& shape) {
        return ArrayND<ScalarType, Dim>::toSize(shape);;
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueTensor<Derived>::grads_impl() const noexcept {
        return GradTensor<Derived, GradOrder>(Base::getDerived());
    }
}
