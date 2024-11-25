/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/GridImpl/RValueGrid.h"

namespace Physica::Core {
    template<class GridType> class LValueGrid;

    template<Grid T>
    class FlattenGrid<T> : public LValueVector<FlattenGrid<T>> {
        using This = FlattenGrid<T>;

        const T& grid;
    public:
        using Base = LValueVector<FlattenGrid<T>>;
        using typename Base::ScalarType;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    public:
        FlattenGrid(const T& grid_) : grid(grid_) {}
        /* Operators */
        FlattenGrid& operator=(const FlattenGrid& obj);
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return grid.getSize(); }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t index);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t index) const;
    };

    template<Grid T>
    FlattenGrid<T>& FlattenGrid<T>::operator=(const FlattenGrid<T>& obj) {
        Base::getDerived() = obj.getDerived();
        return *this;
    }

    template<Grid T>
    __host__ __device__ FlattenGrid<T>::PtrTy FlattenGrid<T>::data_ptr(size_t index) {
        return const_cast<ScalarType*>(const_cast<const This&>(*this).data_ptr(index));
    }

    template<Grid T>
    __host__ __device__ FlattenGrid<T>::ConstPtrTy FlattenGrid<T>::data_ptr(size_t index) const {
        const size_t temp = index / grid.getDimZ();
        const size_t x = temp / grid.getDimY();
        const size_t y = temp % grid.getDimY();
        const size_t z = index % grid.getDimZ();
        return grid.data_ptr({x, y, z});
    }
}

namespace Physica {
    template<Grid T>
    class Traits<FlattenGrid<T>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static size_t SizeAtCompile = Dynamic;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
