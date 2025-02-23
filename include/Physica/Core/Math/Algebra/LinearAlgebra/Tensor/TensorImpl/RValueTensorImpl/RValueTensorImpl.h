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

namespace Physica {
    template<class Derived>
    template<Tensor X>
    void RValueTensor<Derived>::assign(X& x) const {
        forIndexInTensor(getDim(), [this, &x](Index3D index) {
            x(index) = calc(index);
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
    template<Scalar T, bool IsUnitLattice, class Functor>
    inline void RValueTensor<Derived>::forPointInTensor(
            const RValueTensor& grid, const PeriodicCell<T, 3>::LatticeMatrix& lattice, Functor func) {
        return forPointInTensor<T, IsUnitLattice, Functor>(grid.getDim(), lattice, func);
    }

    template<class Derived>
    template<Scalar T, bool IsUnitLattice, class Functor>
    inline void RValueTensor<Derived>::forPointIndexInTensor(
            const RValueTensor& grid, const PeriodicCell<T, 3>::LatticeMatrix& lattice, Functor func) {
        forPointIndexInTensor<T, IsUnitLattice, Functor>(grid.getDim(), lattice, func);
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueTensor<Derived>::grads_impl() const noexcept {
        return GradTensor<Derived, GradOrder>(Base::getDerived());
    }
}
