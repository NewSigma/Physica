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
#include "LGridBlock.h"

namespace Physica::Core {
    template<class GridType> class FlattenGrid;
    /**
     * Notes:
     * Right is positive direction of x.
     * Back is positive direction of y.
     * Top is positive direction of z.
     */
    template<class Derived>
    class LValueGrid : public Utils::CRTPBase<Derived>, public GridBase {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
    public:
        template<class OtherDerived>
        Derived& operator=(const LValueGrid<OtherDerived>& other);
        Derived& operator=(const ScalarType& s);
        [[nodiscard]] ScalarType& operator()(size_t x, size_t y, size_t z) { return Base::getDerived()(x, y, z); }
        [[nodiscard]] const ScalarType& operator()(size_t x, size_t y, size_t z) const { return Base::getDerived()(x, y, z); }
        [[nodiscard]] ScalarType& operator()(Index3D index) { return operator()(index[0], index[1], index[2]); }
        [[nodiscard]] const ScalarType& operator()(Index3D index) const { return operator()(index[0], index[1], index[2]); }
        /* Operations */
        [[nodiscard]] inline LGridBlock<Derived> leftFrontBottomCorner(Index3D cornerIndex);
        [[nodiscard]] inline const LGridBlock<Derived> leftFrontBottomCorner(Index3D cornerIndex) const;
        [[nodiscard]] inline LGridBlock<Derived> rightFrontBottomCorner(Index3D cornerIndex);
        [[nodiscard]] inline const LGridBlock<Derived> rightFrontBottomCorner(Index3D cornerIndex) const;
        [[nodiscard]] inline LGridBlock<Derived> leftBackBottomCorner(Index3D cornerIndex);
        [[nodiscard]] inline const LGridBlock<Derived> leftBackBottomCorner(Index3D cornerIndex) const;
        [[nodiscard]] inline LGridBlock<Derived> rightBackBottomCorner(Index3D cornerIndex);
        [[nodiscard]] inline const LGridBlock<Derived> rightBackBottomCorner(Index3D cornerIndex) const;

        [[nodiscard]] inline LGridBlock<Derived> leftFrontTopCorner(Index3D cornerIndex);
        [[nodiscard]] inline const LGridBlock<Derived> leftFrontTopCorner(Index3D cornerIndex) const;
        [[nodiscard]] inline LGridBlock<Derived> rightFrontTopCorner(Index3D cornerIndex);
        [[nodiscard]] inline const LGridBlock<Derived> rightFrontTopCorner(Index3D cornerIndex) const;
        [[nodiscard]] inline LGridBlock<Derived> leftBackTopCorner(Index3D cornerIndex);
        [[nodiscard]] inline const LGridBlock<Derived> leftBackTopCorner(Index3D cornerIndex) const;
        [[nodiscard]] inline LGridBlock<Derived> rightBackTopCorner(Index3D cornerIndex);
        [[nodiscard]] inline const LGridBlock<Derived> rightBackTopCorner(Index3D cornerIndex) const;

        [[nodiscard]] inline LGridBlock<Derived> block(Index3D from, Index3D count);
        [[nodiscard]] inline const LGridBlock<Derived> block(Index3D from, Index3D count) const;

        [[nodiscard]] FlattenGrid<Derived> flatten() const { return FlattenGrid<Derived>(Base::getDerived()); }
        void resize(Index3D size) { Base::getDerived().resize(size); }
        template<class RandomGenerator> void random_uniform(RandomGenerator& gen);
        template<class RandomGenerator> void random_normal(RandomGenerator& gen);
        /* Getters */
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
        inline static void forPointInGrid(const LValueGrid& grid, const LatticeMatrix& lattice, Functor func);
        template<bool IsUnitLattice, class Functor>
        inline static void forPointIndexInGrid(const LValueGrid& grid, const LatticeMatrix& lattice, Functor func);
    };
}

#include "LValueGridImpl.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/FlattenGrid.h"
