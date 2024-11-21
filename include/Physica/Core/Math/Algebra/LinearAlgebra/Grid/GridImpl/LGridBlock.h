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

namespace Physica::Core {
    template<class Derived> class LValueGrid;
    template<class GridType> class LGridBlock;

    template<Grid T>
    class LGridBlock<T> : public LValueGrid<LGridBlock<T>> {
        using This = LGridBlock<T>;
        using Base = LValueGrid<This>;
        using Index3D = GridBase::Index3D;
    public:
        using typename Base::ScalarType;
    private:
        T& grid;
        Index3D from;
        Index3D count;
    public:
        LGridBlock(T& grid_, Index3D from_, Index3D count_);
        LGridBlock(const LGridBlock&) = delete;
        LGridBlock(LGridBlock&&) noexcept = delete;
        ~LGridBlock() = default;
        /* Operators */
        using Base::operator=;
        LGridBlock& operator=(const LGridBlock& b) { Base::operator=(static_cast<const Base::Base&>(b)); return *this; }
        LGridBlock& operator=(LGridBlock&& b) noexcept { Base::operator=(static_cast<const Base::Base&>(b)); return *this; }
        /* Operations */
        void resize([[maybe_unused]] Index3D size) { assert(size == count && "[Error]: Resize part of a grid is not allowed"); }
        /* Getters */
        [[nodiscard]] size_t getDimX() const noexcept { return count[0]; }
        [[nodiscard]] size_t getDimY() const noexcept { return count[1]; }
        [[nodiscard]] size_t getDimZ() const noexcept { return count[2]; }
        [[nodiscard]] __host__ __device__ inline ScalarType* data_ptr(Index3D index);
        [[nodiscard]] __host__ __device__ inline const ScalarType* data_ptr(Index3D index) const;
    };

    template<Grid T>
    LGridBlock<T>::LGridBlock(T& grid_, Index3D from_, Index3D count_)
            : grid(grid_)
            , from(from_)
            , count(count_) {
        for (int i = 0; i < 3; ++i) {
            assert(from[i] < grid.getDim()[i]);
            assert(from[i] + count[i] <= grid.getDim()[i]);
        }
    }

    template<Grid T>
    __host__ __device__ inline LGridBlock<T>::ScalarType* LGridBlock<T>::data_ptr(Index3D index) {
        return grid.data_ptr({from[0] + index[0], from[1] + index[1], from[2] + index[2]});
    }

    template<Grid T>
    __host__ __device__ inline const LGridBlock<T>::ScalarType* LGridBlock<T>::data_ptr(Index3D index) const {
        return const_cast<This&>(*this).data_ptr(index);
    }
}

namespace Physica {
    template<Grid T>
    class Traits<Core::LGridBlock<T>> {
    public:
        using ScalarType = T::ScalarType;
    };
}
