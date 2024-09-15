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

#include "GridBase.h"

namespace Physica::Core {
    template<class T>
    class GridStorage : public GridBase {
        using This = GridStorage<T>;
        using Base = GridBase;
    public:
        using Index3D = Array<size_t, 3>;
        using ArrayType = Array<T>;
    private:
        ArrayType values;
        size_t dimX;
        size_t dimY;
        size_t dimZ;
    public:
        GridStorage() = default;
        template<class... Args>
        GridStorage(Index3D index, Args&&... args);
        GridStorage(const GridStorage&) = default;
        GridStorage(GridStorage&&) noexcept = default;
        ~GridStorage() = default;
        /* Operators */
        GridStorage& operator=(GridStorage obj) noexcept;
        [[nodiscard]] inline T& operator()(size_t x, size_t y, size_t z) { return *data_ptr({x, y, z}); }
        [[nodiscard]] inline const T& operator()(size_t x, size_t y, size_t z) const { return *data_ptr({x, y, z}); }
        [[nodiscard]] T& operator()(Index3D index) { return *data_ptr(index); }
        [[nodiscard]] const T& operator()(Index3D index) const { return *data_ptr(index); }
        /* Operations */
        template<class Functor> void forIndexInGrid(Functor func) const { Base::forIndexInGrid(getDim(), func); }
        template<class... Args> void resize(Index3D index, Args&&... args);
        void swap(GridStorage& __restrict obj) noexcept;
        /* Iterator */
        auto begin() noexcept { return values.begin(); }
        auto begin() const noexcept { return cbegin(); }
        auto end() noexcept { return values.end(); }
        auto end() const noexcept { return cend(); }
        auto cbegin() const noexcept { return values.cbegin(); }
        auto cend() const noexcept { return values.cend(); }
        auto rbegin() noexcept { return values.rbegin(); }
        auto rend() noexcept { return values.rend(); }
        auto crbegin() const noexcept { return values.crbegin(); }
        auto crend() const noexcept { return values.crend(); }
        /* Getters */
        [[nodiscard]] ArrayType& asArray() noexcept { return values; }
        [[nodiscard]] const ArrayType& asArray() const noexcept { return values; }
        [[nodiscard]] size_t getDimX() const noexcept { return dimX; }
        [[nodiscard]] size_t getDimY() const noexcept { return dimY; }
        [[nodiscard]] size_t getDimZ() const noexcept { return dimZ; }
        [[nodiscard]] Index3D getDim() const noexcept { return {getDimX(), getDimY(), getDimZ()}; }
        [[nodiscard]] size_t getSize() const noexcept { return values.getLength(); }
        [[nodiscard]] __host__ __device__ inline T* data_ptr(Index3D index);
        [[nodiscard]] __host__ __device__ inline const T* data_ptr(Index3D index) const;
    };

    template<class T>
    GridStorage<T>& GridStorage<T>::operator=(GridStorage<T> obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class T>
    template<class... Args>
    GridStorage<T>::GridStorage(Index3D index, Args&&... args)
            : dimX(index[0])
            , dimY(index[1])
            , dimZ(index[2]) {
        values.resize(dimX * dimY * dimZ, std::forward<Args>(args)...);
    }

    template<class T>
    template<class... Args>
    void GridStorage<T>::resize(Index3D index, Args&&... args) {
        dimX = index[0];
        dimY = index[1];
        dimZ = index[2];
        values.resize(dimX * dimY * dimZ, std::forward<Args>(args)...);
    }

    template<class T>
    void GridStorage<T>::swap(GridStorage& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        values.swap(obj.values);
        std::swap(dimX, obj.dimX);
        std::swap(dimY, obj.dimY);
        std::swap(dimZ, obj.dimZ);
    }

    template<class T>
    __host__ __device__ inline T* GridStorage<T>::data_ptr(Index3D index) {
        assert(index[0] < dimX);
        assert(index[1] < dimY);
        assert(index[2] < dimZ);
        return values.data() + (index[0] * dimY + index[1]) * dimZ + index[2];
    }

    template<class T>
    __host__ __device__ inline const T* GridStorage<T>::data_ptr(Index3D index) const {
        return const_cast<This&>(*this).data_ptr(index);
    }
}
