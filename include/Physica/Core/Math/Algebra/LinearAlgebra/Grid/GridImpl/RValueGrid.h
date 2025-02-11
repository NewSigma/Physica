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

#include "Physica/CRTPBase.h"
#include "GridBase.h"

namespace Physica {
    template<class Derived> class RValueGrid;
    template<class Derived> class LValueGrid;
    template<class GridType> class FlattenGrid;
    template<class T> class RealGrid;
    template<class T> class ImagGrid;
    template<class T> class SquaredNormGrid;
    template<class T> class NormGrid;
    template<class T> class ValueGrid;
    template<class T, int GradOrder> class GradGrid;

    template<class T>
    concept Grid = std::derived_from<T, RValueGrid<T>>;

    template<class T>
    concept LGrid = std::derived_from<T, LValueGrid<T>>;

    template<class Derived>
    class RValueGrid : public CRTPBase<RValueGrid<Derived>>, public GridBase {
        using This = RValueGrid<Derived>;
        using Base = CRTPBase<This>;
    public:
        using ScalarType = Traits<Derived>::ScalarType;
        constexpr static bool isForwardDiff = ScalarType::isForwardDiff;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isComplex = ScalarType::isComplex;
    public:
        /* Operations */
        template<Grid T>
        void assign(LValueGrid<T>& other) const;
        template<class Functor> void forIndexInGrid(Functor func) const { GridBase::forIndexInGrid(getDim(), func); }

        auto reals() const noexcept;
        auto imags() const noexcept;
        auto squaredNorms() const noexcept;
        auto norms() const noexcept;
        auto values() const noexcept;
        template<int GradOrder = 1>
        auto grads() const noexcept;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t x, size_t y, size_t z) const { return calc({x, y, z}); }
        [[nodiscard]] ScalarType calc(Index3D index) const { return Base::getDerived().calc(index); }
        [[nodiscard]] size_t getDimX() const noexcept { return Base::getDerived().getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return Base::getDerived().getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return Base::getDerived().getDimZ(); }
        [[nodiscard]] Index3D getDim() const noexcept { return {getDimX(), getDimY(), getDimZ()}; }
        [[nodiscard]] size_t getSize() const noexcept { return getDimX() * getDimY() * getDimZ(); }
        /* Static members */
        using GridBase::forPointInGrid;
        using GridBase::forPointIndexInGrid;
        using GridBase::forIndexInGrid;
        template<Scalar T, bool IsUnitLattice, class Functor>
        inline static void forPointInGrid(
                const RValueGrid& grid, const PeriodicCell<T, 3>::LatticeMatrix& lattice, Functor func);
        template<Scalar T, bool IsUnitLattice, class Functor>
        inline static void forPointIndexInGrid(
                const RValueGrid& grid, const PeriodicCell<T, 3>::LatticeMatrix& lattice, Functor func);
    protected:
        template<int GradOrder>
        auto grads_impl() const noexcept;
    };

    template<class Derived>
    template<Grid T>
    void RValueGrid<Derived>::assign(LValueGrid<T>& other) const {
        forIndexInGrid(getDim(), [this, &other](Index3D index) {
            other(index) = calc(index);
        });
    }

    template<class Derived>
    auto RValueGrid<Derived>::reals() const noexcept {
        return RealGrid<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueGrid<Derived>::imags() const noexcept {
        return ImagGrid<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueGrid<Derived>::squaredNorms() const noexcept {
        return SquaredNormGrid<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueGrid<Derived>::norms() const noexcept {
        return NormGrid<Derived>(Base::getDerived());
    }

    template<class Derived>
    auto RValueGrid<Derived>::values() const noexcept {
        return ValueGrid<Derived>(Base::getDerived());
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueGrid<Derived>::grads() const noexcept {
        if constexpr (isReverseDiff)
            return Base::getDerived().template grads<GradOrder>();
        else
            return grads_impl<GradOrder>();
    }

    template<class Derived>
    template<Scalar T, bool IsUnitLattice, class Functor>
    inline void RValueGrid<Derived>::forPointInGrid(
            const RValueGrid& grid, const PeriodicCell<T, 3>::LatticeMatrix& lattice, Functor func) {
        return forPointInGrid<T, IsUnitLattice, Functor>(grid.getDim(), lattice, func);
    }

    template<class Derived>
    template<Scalar T, bool IsUnitLattice, class Functor>
    inline void RValueGrid<Derived>::forPointIndexInGrid(
            const RValueGrid& grid, const PeriodicCell<T, 3>::LatticeMatrix& lattice, Functor func) {
        forPointIndexInGrid<T, IsUnitLattice, Functor>(grid.getDim(), lattice, func);
    }

    template<class Derived>
    template<int GradOrder>
    auto RValueGrid<Derived>::grads_impl() const noexcept {
        return GradGrid<Derived, GradOrder>(Base::getDerived());
    }
}

namespace Physica {
    template<class T>
    class Traits<RValueGrid<T>> {
    public:
        using Derived = T;
    };
}

#include "GridExpression.h"
#include "GridConvert.h"
