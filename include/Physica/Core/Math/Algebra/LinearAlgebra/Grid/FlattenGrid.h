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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"

namespace Physica::Core {
    template<class GridType> class FlattenGrid;

    namespace Internal {
        template<class T> class Traits;

        template<class GridType>
        class Traits<FlattenGrid<GridType>> {
        public:
            using ScalarType = typename GridType::ScalarType;
            constexpr static size_t SizeAtCompile = Dynamic;
            constexpr static size_t MaxSizeAtCompile = Dynamic;
        };
    }

    template<class GridType>
    class FlattenGrid : public LValueVector<FlattenGrid<GridType>> {
        using This = FlattenGrid<GridType>;

        const GridType& grid;
    public:
        using Base = LValueVector<FlattenGrid<GridType>>;
        using typename Base::ScalarType;
    public:
        FlattenGrid(const LValueGrid<GridType>& grid_) : grid(grid_.getDerived()) {}
        /* Operators */
        FlattenGrid& operator=(const FlattenGrid& obj);
        using Base::operator=;
        [[nodiscard]] ScalarType& operator[](size_t index);
        [[nodiscard]] const ScalarType& operator[](size_t index) const;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return grid.getSize(); }
    };

    template<class GridType>
    FlattenGrid<GridType>& FlattenGrid<GridType>::operator=(const FlattenGrid<GridType>& obj) {
        Base::getDerived() = obj.getDerived();
        return *this;
    }

    template<class GridType>
    typename FlattenGrid<GridType>::ScalarType& FlattenGrid<GridType>::operator[](size_t index) {
        return const_cast<ScalarType&>(const_cast<const This&>(*this).operator[](index));
    }

    template<class GridType>
    const typename FlattenGrid<GridType>::ScalarType& FlattenGrid<GridType>::operator[](size_t index) const {
        const size_t temp = index / grid.getDimZ();
        const size_t x = temp / grid.getDimY();
        const size_t y = temp % grid.getDimY();
        const size_t z = index % grid.getDimZ();
        return grid(x, y, z);
    }
}
