/*
 * Copyright 2023 Weibo He.
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

#include "RValueGrid.h"
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
    class LValueGrid : public RValueGrid<Derived> {
        using Base = RValueGrid<Derived>;
    public:
        using typename Base::ScalarType;
        using typename Base::Index3D;
    public:
        template<class OtherDerived>
        Derived& operator=(const RValueGrid<OtherDerived>& other);
        Derived& operator=(const ScalarType& s);
        [[nodiscard]] ScalarType& operator()(size_t x, size_t y, size_t z) { return *data_ptr({x, y, z}); }
        [[nodiscard]] const ScalarType& operator()(size_t x, size_t y, size_t z) const { return *data_ptr({x, y, z}); }
        [[nodiscard]] ScalarType& operator()(Index3D index) { return *data_ptr(index); }
        [[nodiscard]] const ScalarType& operator()(Index3D index) const { return *data_ptr(index); }
        template<class T> void operator+=(const ScalarBase<T>& s) { (*this) = (*this) + s.getDerived(); }
        template<class T> void operator-=(const ScalarBase<T>& s) { (*this) = (*this) - s.getDerived(); }
        template<class T> void operator*=(const ScalarBase<T>& s) { (*this) = (*this) * s.getDerived(); }
        template<class T> void operator/=(const ScalarBase<T>& s) { (*this) = (*this) / s.getDerived(); }
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
        [[nodiscard]] ScalarType calc(Index3D index) const { return operator()(index); }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(Index3D index) { return Base::getDerived().data_ptr(index); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(Index3D index) const { return Base::getDerived().data_ptr(index); }
        using Base::getDimX;
        using Base::getDimY;
        using Base::getDimZ;
        using Base::getDim;
        /* Static members */
        using Base::forPointInGrid;
        using Base::forPointIndexInGrid;
        using Base::forIndexInGrid;
    };

    template<class Derived, class OtherDerived>
    inline void operator+=(LValueGrid<Derived>& g1, const RValueGrid<OtherDerived>& g2) {
        g1.getDerived() = g1.getDerived() + g2.getDerived();
    }

    template<class Derived, class OtherDerived>
    inline void operator-=(LValueGrid<Derived>& g1, const RValueGrid<OtherDerived>& g2) {
        g1.getDerived() = g1.getDerived() - g2.getDerived();
    }
}

#include "LValueGridImpl.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/FlattenGrid.h"
