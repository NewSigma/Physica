/*
 * Copyright 2021 WeiBo He.
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
    template<class ScalarType>
    class Grid3D {
    public:
        using LatticeMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Element, 3, 3>;
        using Dim = std::tuple<size_t, size_t, size_t>;
    private:
        LatticeMatrix lattice;
        Vector<ScalarType> values;
        size_t dimX;
        size_t dimY;
        size_t dimZ;
    public:
        template<class MatrixType>
        Grid3D(const RValueMatrix<MatrixType>& lattice_, size_t dimX_, size_t dimY_, size_t dimZ_);
        /* Operators */
        [[nodiscard]] ScalarType& operator[](size_t index);
        [[nodiscard]] const ScalarType& operator[](size_t index) const;
        [[nodiscard]] ScalarType& operator()(size_t x, size_t y, size_t z);
        [[nodiscard]] const ScalarType& operator()(size_t x, size_t y, size_t z) const;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] ScalarType getVolume() const noexcept;
        [[nodiscard]] ScalarType getUnitVolume() const noexcept;
        [[nodiscard]] Vector<ScalarType>& asVector() noexcept { return values; }
        [[nodiscard]] const Vector<ScalarType>& asVector() const noexcept { return values; }
        [[nodiscard]] size_t getSize() const noexcept;
        [[nodiscard]] Dim getDim() const noexcept;
        [[nodiscard]] Vector<ScalarType, 3> dimToPos(Dim dim) const;
        [[nodiscard]] Vector<ScalarType, 3> indexToPos(size_t index) const;
        [[nodiscard]] size_t dimToIndex(size_t x, size_t y, size_t z) const noexcept;
        [[nodiscard]] Dim indexToDim(size_t index) const noexcept;
    };

    template<class ScalarType>
    template<class MatrixType>
    Grid3D<ScalarType>::Grid3D(const RValueMatrix<MatrixType>& lattice_, size_t dimX_, size_t dimY_, size_t dimZ_)
            : lattice(lattice_)
            , values(dimX_ * dimY_ * dimZ_)
            , dimX(dimX_)
            , dimY(dimY_)
            , dimZ(dimZ_) {}

    template<class ScalarType>
    ScalarType& Grid3D<ScalarType>::operator[](size_t index) {
        return values[index];
    }

    template<class ScalarType>
    const ScalarType& Grid3D<ScalarType>::operator[](size_t index) const {
        return values[index];
    }

    template<class ScalarType>
    ScalarType& Grid3D<ScalarType>::operator()(size_t x, size_t y, size_t z) {
        return values[dimToIndex(x, y, z)];
    }

    template<class ScalarType>
    const ScalarType& Grid3D<ScalarType>::operator()(size_t x, size_t y, size_t z) const {
        return values[dimToIndex(x, y, z)];
    }

    template<class ScalarType>
    ScalarType Grid3D<ScalarType>::getVolume() const noexcept {
        return abs((lattice.row(0).crossProduct(lattice.row(1))).compute() * lattice.row(2));
    }

    template<class ScalarType>
    ScalarType Grid3D<ScalarType>::getUnitVolume() const noexcept {
        return getVolume() / ScalarType((dimX - 1) * (dimY - 1) * (dimZ - 1));
    }

    template<class ScalarType>
    size_t Grid3D<ScalarType>::getSize() const noexcept {
        return dimX * dimY * dimZ;
    }

    template<class ScalarType>
    typename Grid3D<ScalarType>::Dim Grid3D<ScalarType>::getDim() const noexcept {
        return std::make_tuple(dimX, dimY, dimZ);
    }

    template<class ScalarType>
    Vector<ScalarType, 3> Grid3D<ScalarType>::dimToPos(Dim dim) const {
        auto[x, y, z] = dim;
        const ScalarType factor_x = ScalarType(x) / ScalarType(dimX - 1);
        const ScalarType factor_y = ScalarType(y) / ScalarType(dimY - 1);
        const ScalarType factor_z = ScalarType(z) / ScalarType(dimZ - 1);
        return lattice.row(0) * factor_x + lattice.row(1) * factor_y + lattice(2) * factor_z;
    }

    template<class ScalarType>
    Vector<ScalarType, 3> Grid3D<ScalarType>::indexToPos(size_t index) const {
        return dimToPos(indexToDim(index));
    }

    template<class ScalarType>
    size_t Grid3D<ScalarType>::dimToIndex(size_t x, size_t y, size_t z) const noexcept {
        const size_t dimXY = dimX * dimY;
        return z * dimXY + y * dimX + x;
    }

    template<class ScalarType>
    typename Grid3D<ScalarType>::Dim Grid3D<ScalarType>::indexToDim(size_t index) const noexcept {
        const size_t dimXY = dimX * dimY;

        const size_t z = index / dimXY;
        index -= dimXY * z;
        const size_t y = index / dimX;
        index -= dimX * y;
        const size_t x = index;

        return {x, y, z};
    }
}
