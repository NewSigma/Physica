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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/CrossProduct.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/PhyConst.h"

namespace Physica::Core {
    /**
     * \tparam isSigned
     * If true, the grid is [-dimX, dimX] * [-dimY, dimY] * [-dimZ, dimZ]
     * If false, the grid is [0, dimX - 1] * [0, dimY - 1] * [0, dimZ - 1]
     */
    template<class ScalarType, bool isSigned>
    class Grid3D {
        using RealType = typename ScalarType::RealType;
    public:
        using IntType = typename std::conditional<isSigned, ssize_t, size_t>::type;
        using LatticeMatrix = typename CrystalCell::LatticeMatrix;
        using Dim = std::tuple<IntType, IntType, IntType>;
    private:
        LatticeMatrix lattice;
        Vector<ScalarType> values;
        size_t dimX;
        size_t dimY;
        size_t dimZ;
    public:
        Grid3D() = default;
        template<class MatrixType>
        Grid3D(const RValueMatrix<MatrixType>& lattice_, size_t dimX_, size_t dimY_, size_t dimZ_);
        Grid3D(const Grid3D&) = default;
        Grid3D(Grid3D&&) noexcept = default;
        ~Grid3D() = default;
        /* Operators */
        Grid3D& operator=(Grid3D grid);
        [[nodiscard]] ScalarType& operator[](size_t index);
        [[nodiscard]] const ScalarType& operator[](size_t index) const;
        [[nodiscard]] ScalarType& operator()(IntType x, IntType y, IntType z);
        [[nodiscard]] const ScalarType& operator()(IntType x, IntType y, IntType z) const;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] size_t getDimX() const noexcept { return dimX; }
        [[nodiscard]] size_t getDimY() const noexcept { return dimY; }
        [[nodiscard]] size_t getDimZ() const noexcept { return dimZ; }
        [[nodiscard]] ScalarType getVolume() const noexcept;
        [[nodiscard]] ScalarType getUnitVolume() const noexcept;
        [[nodiscard]] Vector<ScalarType>& asVector() noexcept { return values; }
        [[nodiscard]] const Vector<ScalarType>& asVector() const noexcept { return values; }
        [[nodiscard]] size_t getSize() const noexcept;
        [[nodiscard]] Dim getDim() const noexcept;
        [[nodiscard]] Vector<RealType, 3> dimToPos(Dim dim) const;
        [[nodiscard]] Vector<RealType, 3> indexToPos(size_t index) const;
        [[nodiscard]] size_t dimToIndex(IntType x, IntType y, IntType z) const noexcept;
        [[nodiscard]] Dim indexToDim(size_t index) const noexcept;
        /* Helpers */
        void swap(Grid3D& grid) noexcept;
        [[nodiscard]] static Grid3D gridFromCutEnergy(RealType cutEnergy, LatticeMatrix reciprocalLattice);
    };

    template<class ScalarType, bool isSigned>
    template<class MatrixType>
    Grid3D<ScalarType, isSigned>::Grid3D(const RValueMatrix<MatrixType>& lattice_, size_t dimX_, size_t dimY_, size_t dimZ_)
            : lattice(lattice_)
            , dimX(dimX_)
            , dimY(dimY_)
            , dimZ(dimZ_) {
        values.resize(getSize());
    }

    template<class ScalarType, bool isSigned>
    Grid3D<ScalarType, isSigned>& Grid3D<ScalarType, isSigned>::operator=(Grid3D grid) {
        swap(grid);
        return *this;
    }

    template<class ScalarType, bool isSigned>
    ScalarType& Grid3D<ScalarType, isSigned>::operator[](size_t index) {
        return values[index];
    }

    template<class ScalarType, bool isSigned>
    const ScalarType& Grid3D<ScalarType, isSigned>::operator[](size_t index) const {
        return values[index];
    }

    template<class ScalarType, bool isSigned>
    ScalarType& Grid3D<ScalarType, isSigned>::operator()(IntType x, IntType y, IntType z) {
        return values[dimToIndex(x, y, z)];
    }

    template<class ScalarType, bool isSigned>
    const ScalarType& Grid3D<ScalarType, isSigned>::operator()(IntType x, IntType y, IntType z) const {
        return values[dimToIndex(x, y, z)];
    }

    template<class ScalarType, bool isSigned>
    ScalarType Grid3D<ScalarType, isSigned>::getVolume() const noexcept {
        return abs((lattice.row(0).crossProduct(lattice.row(1))).compute() * lattice.row(2));
    }

    template<class ScalarType, bool isSigned>
    ScalarType Grid3D<ScalarType, isSigned>::getUnitVolume() const noexcept {
        if constexpr (isSigned)
            return getVolume() / ScalarType((dimX - 1) * (dimY - 1) * (dimZ - 1) * 8);
        else
            return getVolume() / ScalarType((dimX - 1) * (dimY - 1) * (dimZ - 1));
    }

    template<class ScalarType, bool isSigned>
    size_t Grid3D<ScalarType, isSigned>::getSize() const noexcept {
        if constexpr (isSigned)
            return (2 * dimX + 1) * (2 * dimY + 1) * (2 * dimZ + 1);
        else
            return dimX * dimY * dimZ;
    }

    template<class ScalarType, bool isSigned>
    typename Grid3D<ScalarType, isSigned>::Dim Grid3D<ScalarType, isSigned>::getDim() const noexcept {
        return std::make_tuple(dimX, dimY, dimZ);
    }

    template<class ScalarType, bool isSigned>
    Vector<typename Grid3D<ScalarType, isSigned>::RealType, 3> Grid3D<ScalarType, isSigned>::dimToPos(Dim dim) const {
        auto[x, y, z] = dim;
        if constexpr (isSigned) {
            return lattice.row(0).asVector() * RealType(x) +
                   lattice.row(1).asVector() * RealType(y) +
                   lattice.row(2).asVector() * RealType(z);
        }
        else {
            const RealType factor_x = RealType(x) / RealType(dimX - 1);
            const RealType factor_y = RealType(y) / RealType(dimY - 1);
            const RealType factor_z = RealType(z) / RealType(dimZ - 1);
            return lattice.row(0).asVector() * factor_x + lattice.row(1).asVector() * factor_y + lattice.row(2).asVector() * factor_z;
        }
    }

    template<class ScalarType, bool isSigned>
    Vector<typename Grid3D<ScalarType, isSigned>::RealType, 3> Grid3D<ScalarType, isSigned>::indexToPos(size_t index) const {
        return dimToPos(indexToDim(index));
    }

    template<class ScalarType, bool isSigned>
    size_t Grid3D<ScalarType, isSigned>::dimToIndex(IntType x, IntType y, IntType z) const noexcept {
        if constexpr (isSigned) {
            const size_t sizeX = 2 * dimX + 1;
            const size_t sizeXY = sizeX * (2 * dimY + 1);
            
            const size_t normalized_x = x + static_cast<ssize_t>(dimX);
            const size_t normalized_y = y + static_cast<ssize_t>(dimY);
            const size_t normalized_z = z + static_cast<ssize_t>(dimZ);
            return normalized_z * sizeXY + normalized_y * sizeX + normalized_x;
        }
        else {
            const size_t dimXY = dimX * dimY;
            return z * dimXY + y * dimX + x;
        }
    }

    template<class ScalarType, bool isSigned>
    typename Grid3D<ScalarType, isSigned>::Dim Grid3D<ScalarType, isSigned>::indexToDim(size_t index) const noexcept {
        if constexpr (isSigned) {
            const size_t sizeX = 2 * dimX + 1;
            const size_t sizeXY = sizeX * (2 * dimY + 1);

            const size_t normalized_z = index / sizeXY;
            index -= sizeXY * normalized_z;
            const size_t normalized_y = index / sizeX;
            index -= sizeX * normalized_y;
            const size_t normalized_x = index;

            return {static_cast<ssize_t>(normalized_x) - static_cast<ssize_t>(dimX),
                    static_cast<ssize_t>(normalized_y) - static_cast<ssize_t>(dimY),
                    static_cast<ssize_t>(normalized_z) - static_cast<ssize_t>(dimZ)};
        }
        else {
            const size_t dimXY = dimX * dimY;

            const size_t z = index / dimXY;
            index -= dimXY * z;
            const size_t y = index / dimX;
            index -= dimX * y;
            const size_t x = index;

            return {x, y, z};
        }
    }

    template<class ScalarType, bool isSigned>
    void Grid3D<ScalarType, isSigned>::swap(Grid3D& grid) noexcept {
        std::swap(lattice, grid.lattice);
        std::swap(values, grid.values);
        std::swap(dimX, grid.dimX);
        std::swap(dimY, grid.dimY);
        std::swap(dimZ, grid.dimZ);
    }

    template<class ScalarType, bool isSigned>
    Grid3D<ScalarType, isSigned> Grid3D<ScalarType, isSigned>::gridFromCutEnergy(RealType cutEnergy, LatticeMatrix reciprocalLattice) {
        const auto factor = RealType(2 * PhyConst<AU>::electronMass / PhyConst<AU>::reducedPlanck / PhyConst<AU>::reducedPlanck);
        const RealType maxMoment = sqrt(factor * cutEnergy);
        const auto dimX = size_t((maxMoment / reciprocalLattice.row(0).norm()).getTrivial());
        const auto dimY = size_t((maxMoment / reciprocalLattice.row(1).norm()).getTrivial());
        const auto dimZ = size_t((maxMoment / reciprocalLattice.row(2).norm()).getTrivial());
        return Grid3D<ScalarType, isSigned>(reciprocalLattice, dimX, dimY, dimZ);
    }
}

namespace std {
    template<class ScalarType, bool isSigned>
    inline void swap(Physica::Core::Grid3D<ScalarType, isSigned>& grid1,
                     Physica::Core::Grid3D<ScalarType, isSigned>& grid2) noexcept {
        grid1.swap(grid2);
    }
}
