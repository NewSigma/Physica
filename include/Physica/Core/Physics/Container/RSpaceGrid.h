/*
 * Copyright 2021-2022 WeiBo He.
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

#include <type_traits>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/CrossProduct.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    template<class T>
    class RSpaceGrid {
        constexpr static bool isScalar = is_scalar<T>::value;
    public:
        using Container = typename std::conditional<isScalar, Vector<T>, Utils::Array<T>>::type;
        using ValueType = T;
        using Index3D = Utils::Array<size_t, 3>;
        using LatticeMatrix = typename CrystalCell::LatticeMatrix;
        using VectorType = Vector<typename LatticeMatrix::ScalarType, 3>;
    private:
        Container values;
        size_t dimX;
        size_t dimY;
        size_t dimZ;
    public:
        RSpaceGrid() = default;
        template<class... Args>
        RSpaceGrid(Index3D index, Args... args);
        RSpaceGrid(const RSpaceGrid&) = default;
        RSpaceGrid(RSpaceGrid&&) noexcept = default;
        ~RSpaceGrid() = default;
        /* Operators */
        RSpaceGrid& operator=(RSpaceGrid grid) noexcept;
        [[nodiscard]] T& operator()(size_t x, size_t y, size_t z);
        [[nodiscard]] const T& operator()(size_t x, size_t y, size_t z) const;
        [[nodiscard]] T& operator()(Index3D index) { return this->operator()(index[0], index[1], index[2]); }
        [[nodiscard]] const T& operator()(Index3D index) const { return this->operator()(index[0], index[1], index[2]); }
        template<class T1>
        friend std::ostream& operator<<(std::ostream& os, const RSpaceGrid<T1>& grid);
        template<class T1>
        friend std::istream& operator>>(std::istream& is, RSpaceGrid<T1>& grid);
        /* Operations */
        template<class... Args>
        void resize(Index3D index, Args... args);
        void swap(RSpaceGrid& grid) noexcept;
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
        [[nodiscard]] Container& flatten() noexcept { return values; }
        [[nodiscard]] const Container& flatten() const noexcept { return values; }
        [[nodiscard]] size_t getSize() const noexcept { return values.getLength(); }
        [[nodiscard]] size_t getDimX() const noexcept { return dimX; }
        [[nodiscard]] size_t getDimY() const noexcept { return dimY; }
        [[nodiscard]] size_t getDimZ() const noexcept { return dimZ; }
        [[nodiscard]] Index3D getDim() const noexcept { return {getDimX(), getDimY(), getDimZ()}; }
        /* Static members */
        template<class Functor>
        static void forPointInGrid(const RSpaceGrid<T>& grid, const LatticeMatrix& lattice, Functor func);
        template<class Functor>
        static void forPointIndexInGrid(const RSpaceGrid<T>& grid, const LatticeMatrix& lattice, Functor func);
        template<class Functor>
        static void forIndexInGrid(Index3D dim, Functor func);
    protected:
        RSpaceGrid(Container values_, size_t dimX_, size_t dimY_, size_t dimZ_);
    };

    template<class T>
    template<class... Args>
    RSpaceGrid<T>::RSpaceGrid(Index3D index, Args... args)
            : dimX(index[0])
            , dimY(index[1])
            , dimZ(index[2]) {
        values.resize(dimX * dimY * dimZ, std::forward<Args>(args)...);
    }

    template<class T>
    RSpaceGrid<T>::RSpaceGrid(Container values_, size_t dimX_, size_t dimY_, size_t dimZ_)
            : values(std::move(values_))
            , dimX(dimX_)
            , dimY(dimY_)
            , dimZ(dimZ_) {
        assert(values.getLength() == dimX * dimY * dimZ);
    }

    template<class T>
    RSpaceGrid<T>& RSpaceGrid<T>::operator=(RSpaceGrid grid) noexcept {
        swap(grid);
        return *this;
    }

    template<class T>
    T& RSpaceGrid<T>::operator()(size_t x, size_t y, size_t z) {
        return values[x * dimY * dimZ + y * dimZ + z];
    }

    template<class T>
    const T& RSpaceGrid<T>::operator()(size_t x, size_t y, size_t z) const {
        return values[x * dimY * dimZ + y * dimZ + z];
    }

    template<class T>
    std::ostream& operator<<(std::ostream& os, const RSpaceGrid<T>& grid) {
        using Index3D = typename RSpaceGrid<T>::Index3D;
        const Index3D dim = grid.getDim();
        os.write(reinterpret_cast<const char*>(&dim), sizeof(Index3D));
        os << grid.flatten();
        return os;
    }

    template<class T>
    std::istream& operator>>(std::istream& is, RSpaceGrid<T>& grid) {
        using Index3D = typename RSpaceGrid<T>::Index3D;
        Index3D dim;
        is.read(reinterpret_cast<char*>(&dim), sizeof(Index3D));
        grid.resize(dim);
        is >> grid.flatten();
        return is;
    }

    template<class T>
    template<class... Args>
    void RSpaceGrid<T>::resize(Index3D index, Args... args) {
        dimX = index[0];
        dimY = index[1];
        dimZ = index[2];
        values.resize(dimX * dimY * dimZ, std::forward<Args>(args)...);
    }

    template<class T>
    void RSpaceGrid<T>::swap(RSpaceGrid& grid) noexcept {
        using std::swap;
        swap(values, grid.values);
        swap(dimX, grid.dimX);
        swap(dimY, grid.dimY);
        swap(dimZ, grid.dimZ);
    }

    template<class T>
    template<class Functor>
    void RSpaceGrid<T>::forPointInGrid(const RSpaceGrid<T>& grid, const LatticeMatrix& lattice, Functor func) {
        using ScalarType = typename VectorType::ScalarType;
        LatticeMatrix sub_lattice{};
        auto a1 = sub_lattice.row(0);
        a1 = lattice.row(0).asVector() * reciprocal(ScalarType(grid.getDimX()));
        auto a2 = sub_lattice.row(1);
        a2 = lattice.row(1).asVector() * reciprocal(ScalarType(grid.getDimY()));
        auto a3 = sub_lattice.row(2);
        a3 = lattice.row(2).asVector() * reciprocal(ScalarType(grid.getDimZ()));

        VectorType v1, v2, v3;
        for (size_t x = 0; x < grid.getDimX(); ++x) {
            v1 = ScalarType(x) * a1.asVector();
            for (size_t y = 0; y < grid.getDimY(); ++y) {
                v2 = v1 + ScalarType(y) * a2.asVector();
                for (size_t z = 0; z < grid.getDimZ(); ++z) {
                    v3 = v2 + ScalarType(z) * a3.asVector();
                    func(v3);
                }
            }
        }
    }

    template<class T>
    template<class Functor>
    void RSpaceGrid<T>::forPointIndexInGrid(const RSpaceGrid<T>& grid, const LatticeMatrix& lattice, Functor func) {
        using ScalarType = typename VectorType::ScalarType;
        LatticeMatrix sub_lattice{};
        auto a1 = sub_lattice.row(0);
        a1 = lattice.row(0).asVector() * reciprocal(ScalarType(grid.getDimX()));
        auto a2 = sub_lattice.row(1);
        a2 = lattice.row(1).asVector() * reciprocal(ScalarType(grid.getDimY()));
        auto a3 = sub_lattice.row(2);
        a3 = lattice.row(2).asVector() * reciprocal(ScalarType(grid.getDimZ()));

        VectorType v1, v2, v3;
        for (size_t x = 0; x < grid.getDimX(); ++x) {
            v1 = ScalarType(x) * a1.asVector();
            for (size_t y = 0; y < grid.getDimY(); ++y) {
                v2 = v1 + ScalarType(y) * a2.asVector();
                for (size_t z = 0; z < grid.getDimZ(); ++z) {
                    v3 = v2 + ScalarType(z) * a3.asVector();
                    func(v3, Index3D{x, y, z});
                }
            }
        }
    }

    template<class T>
    template<class Functor>
    void RSpaceGrid<T>::forIndexInGrid(Index3D dim, Functor func) {
        for (size_t x = 0; x < dim[0]; ++x)
            for (size_t y = 0; y < dim[1]; ++y)
                for (size_t z = 0; z < dim[2]; ++z)
                    func(Index3D{x, y, z});
    }

    template<class T>
    inline void swap(RSpaceGrid<T>& grid1, RSpaceGrid<T>& grid2) noexcept {
        grid1.swap(grid2);
    }
}
