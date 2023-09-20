/*
 * Copyright 2023 WeiBo He.
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

#include "Physica/Utils/Template/CRTPBase.h"
#include "GridBase.h"

namespace Physica::Core {
    template<class Derived> class LValueGrid;

    template<class Derived>
    class RValueGrid : public Utils::CRTPBase<Derived>, public GridBase {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
        using typename GridBase::Index3D;
    public:
        /* Operations */
        template<class OtherDerived>
        void assignTo(LValueGrid<OtherDerived>& other) const;
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
        template<bool IsUnitLattice, class Functor>
        inline static void forPointInGrid(const RValueGrid& grid, const LatticeMatrix& lattice, Functor func);
        template<bool IsUnitLattice, class Functor>
        inline static void forPointIndexInGrid(const RValueGrid& grid, const LatticeMatrix& lattice, Functor func);
    };

    template<class Derived>
    template<class OtherDerived>
    void RValueGrid<Derived>::assignTo(LValueGrid<OtherDerived>& other) const {
        forIndexInGrid(getDim(), [this, &other](Index3D index) {
            other(index) = calc(index);
        });
    }

    template<class Derived>
    template<bool IsUnitLattice, class Functor>
    inline void RValueGrid<Derived>::forPointInGrid(const RValueGrid& grid, const LatticeMatrix& lattice, Functor func) {
        return forPointInGrid<IsUnitLattice, Functor>(grid.getDim(), lattice, func);
    }

    template<class Derived>
    template<bool IsUnitLattice, class Functor>
    inline void RValueGrid<Derived>::forPointIndexInGrid(const RValueGrid& grid, const LatticeMatrix& lattice, Functor func) {
        forPointIndexInGrid<IsUnitLattice, Functor>(grid.getDim(), lattice, func);
    }
}

#include "GridExpression.h"
#include "GridConvert.h"
