/*
 * Copyright 2021-2024 Weibo He.
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

#include "GridImpl/GridStorage.h"
#include "GridImpl/LValueGrid.h"

namespace Physica::Core {
    template<Scalar T>
    class RSpaceGrid : public LValueGrid<RSpaceGrid<T>>, private GridStorage<T> {
        using This = RSpaceGrid<T>;
        using Base = LValueGrid<This>;
        using Storage = GridStorage<T>;
    public:
        using Index3D = typename Storage::Index3D;
    public:
        RSpaceGrid() = default;
        template<class... Args>
        RSpaceGrid(Index3D index, Args&&... args);
        RSpaceGrid(const RSpaceGrid&) = default;
        RSpaceGrid(RSpaceGrid&&) noexcept = default;
        ~RSpaceGrid() = default;
        /* Operators */
        using Base::operator=;
        RSpaceGrid& operator=(RSpaceGrid grid) noexcept;
        using Storage::operator();
        template<class U> friend std::ostream& operator<<(std::ostream& os, const RSpaceGrid<U>& grid);
        template<class U> friend std::istream& operator>>(std::istream& is, RSpaceGrid<U>& grid);
        /* Operations */
        using Base::random_normal;
        template<class... Args>
        inline void resize(Index3D index, Args&&... args);
        inline void swap(RSpaceGrid& __restrict grid) noexcept;
        /* Getters */
        using Storage::data_ptr;
        using Storage::getDim;
        using Storage::getDimX;
        using Storage::getDimY;
        using Storage::getDimZ;
        using Storage::getSize;
        /* Static members */
        using Base::forIndexInGrid;
        using Base::forPointIndexInGrid;
        using Base::forPointInGrid;
        template<class RandomType>
        static RSpaceGrid random_uniform(Index3D size, RandomType& gen);
        template<class RandomType>
        static RSpaceGrid random_normal(Index3D size, RandomType& gen);
    };

    template<Scalar T>
    template<class... Args>
    RSpaceGrid<T>::RSpaceGrid(Index3D index, Args&&... args) : Storage(index, std::forward<Args>(args)...) {}

    template<Scalar T>
    RSpaceGrid<T>& RSpaceGrid<T>::operator=(RSpaceGrid grid) noexcept {
        swap(grid);
        return *this;
    }

    template<Scalar T>
    std::ostream& operator<<(std::ostream& os, const RSpaceGrid<T>& grid) {
        using Index3D = typename RSpaceGrid<T>::Index3D;
        const Index3D dim = grid.getDim();
        os.write(reinterpret_cast<const char*>(&dim), sizeof(Index3D));
        os.write(reinterpret_cast<const char*>(grid.asArray().data()), grid.getSize() * sizeof(T));
        return os;
    }

    template<Scalar T>
    std::istream& operator>>(std::istream& is, RSpaceGrid<T>& grid) {
        using Index3D = typename RSpaceGrid<T>::Index3D;
        Index3D dim;
        is.read(reinterpret_cast<char*>(&dim), sizeof(Index3D));
        grid.resize(dim);
        is.read(reinterpret_cast<char*>(grid.asArray().data()), grid.getSize() * sizeof(T));
        return is;
    }

    template<Scalar T>
    template<class... Args>
    inline void RSpaceGrid<T>::resize(Index3D index, Args&&... args) {
        Storage::resize(index, std::forward<Args>(args)...);
    }

    template<Scalar T>
    inline void RSpaceGrid<T>::swap(RSpaceGrid& __restrict grid) noexcept {
        assert(this != &grid && "[Error]: Self swap is likely a bug");
        Storage::swap(grid);
    }

    template<Scalar T>
    template<class RandomType>
    inline RSpaceGrid<T> RSpaceGrid<T>::random_uniform(Index3D size, RandomType& gen) {
        auto result = RSpaceGrid<T>(size);
        result.flatten().random_uniform(gen);
        return result;
    }

    template<Scalar T>
    template<class RandomType>
    inline RSpaceGrid<T> RSpaceGrid<T>::random_normal(Index3D size, RandomType& gen) {
        auto result = RSpaceGrid<T>(size);
        result.flatten().random_normal(gen);
        return result;
    }

    template<Scalar T>
    inline void swap(RSpaceGrid<T>& __restrict grid1, RSpaceGrid<T>& __restrict grid2) noexcept {
        grid1.swap(grid2);
    }
}

namespace Physica {
    template<Scalar T>
    class Traits<RSpaceGrid<T>> {
    public:
        using ScalarType = T;
    };
}
