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

namespace Physica::Core {
    template<class GridType> class FlattenGrid;

    template<class Derived>
    class LValueGrid : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
    public:
        template<class OtherDerived>
        Derived& operator=(const LValueGrid<OtherDerived>& other);
        [[nodiscard]] ScalarType& operator()(size_t x, size_t y, size_t z) { return Base::getDerived()(x, y, z); }
        [[nodiscard]] const ScalarType& operator()(size_t x, size_t y, size_t z) const { return Base::getDerived()(x, y, z); }
        /* Operations */
        [[nodiscard]] FlattenGrid<Derived> flatten() const { return FlattenGrid<Derived>(Base::getDerived()); }
        /* Getters */
        [[nodiscard]] size_t getDimX() const noexcept { return Base::getDerived().getDimX(); }
        [[nodiscard]] size_t getDimY() const noexcept { return Base::getDerived().getDimY(); }
        [[nodiscard]] size_t getDimZ() const noexcept { return Base::getDerived().getDimZ(); }
        [[nodiscard]] size_t getSize() const noexcept { return getDimX() * getDimY() * getDimZ(); }
    };

    template<class Derived>
    template<class OtherDerived>
    Derived& LValueGrid<Derived>::operator=(const LValueGrid<OtherDerived>& other) {
        assert(getDimX() == other.getDimX());
        assert(getDimY() == other.getDimY());
        assert(getDimZ() == other.getDimZ());
        for (size_t i = 0; i < other.getDimX(); ++i)
            for (size_t j = 0; j < other.getDimY(); ++j)
                for (size_t k = 0; k < other.getDimZ(); ++k)
                    Base::getDerived()(i, j, k) = other.getDerived()(i, j, k);
    }

    template<class Derived>
    void operator*=(LValueGrid<Derived>& grid, typename Derived::ScalarType factor) {
        for (size_t i = 0; i < grid.getDimX(); ++i)
            for (size_t j = 0; j < grid.getDimY(); ++j)
                for (size_t k = 0; k < grid.getDimZ(); ++k)
                    grid(i, j, k) *= factor;
    }
}

#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/FlattenGrid.h"
