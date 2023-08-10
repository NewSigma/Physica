/*
 * Copyright 2022 WeiBo He.
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

#include <forward_list>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/GridImpl/GridStorage.h"
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType>
    class CellList {
    public:
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using CellGrid = GridStorage<std::forward_list<size_t>>;
        using Vector3D = Vector<PosScalarType, 3>;
        using Index3D = typename CellGrid::Index3D;
    private:
        LatticeMatrix lattice;
        PositionMatrix directPos;
        CellGrid cellGrid;
        ScalarType cutoff;
        Utils::Array<Index3D> atomCellMap;
        Utils::Array<Vector3D> neighShifts;
    public:
        CellList(LatticeMatrix lattice_, PositionMatrix pos, ScalarType cutoff_);
        CellList(const MDCellType& mdCell, ScalarType cutoff_);
        CellList(const CellList&) = default;
        CellList(CellList&&) noexcept = default;
        ~CellList() = default;
        /* Operators */
        CellList& operator=(CellList list) noexcept;
        [[nodiscard]] const std::forward_list<size_t>& operator()(Index3D index) const { return cellGrid(index); }
        /* Operations */
        void update(const MDCellType& mdCell);
        void swap(CellList& list) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getCutoff() const noexcept { return cutoff; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return atomCellMap.getLength(); }
        [[nodiscard]] const Utils::Array<Index3D>& getAtomCellMap() const noexcept { return atomCellMap; }
        [[nodiscard]] size_t getNumCell() const noexcept { return cellGrid.getSize(); }
        /* Static members */
        template<class Functor>
        void forCellInList(Functor func) const;
        template<class Functor>
        void forNeighInRange(Index3D centerCell, Functor func) const;
        template<class Functor>
        void forReducedNeighInRange(Index3D centerCell, Functor func) const;
    private:
        Index3D posToIndex(size_t atomId) const;
        /* Helpers */
        void updateNeighShifts();
        template<size_t DimID>
        inline int findNeighbor(size_t centerIndex, int deltaIndex, Index3D& neighborIndex) const;
        /* Static members */
        static Index3D makeGridDim(const LatticeMatrix& lattice, ScalarType cutoff);
    };

    template<class ScalarType, class PosScalarType>
    CellList<ScalarType, PosScalarType>::CellList(LatticeMatrix lattice_, PositionMatrix pos, ScalarType cutoff_)
            : lattice(std::move(lattice_))
            , directPos(std::move(pos))
            , cutoff(cutoff_)
            , neighShifts(3 * 3 * 3) {
        atomCellMap.resize(directPos.getRow());
        const auto dim = makeGridDim(lattice, cutoff);
        cellGrid.resize(dim);

        PeriodicCell<PosScalarType, 3>::toDirect(directPos, lattice);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Index3D index = posToIndex(i);
            atomCellMap[i] = index;
            cellGrid(index).push_front(i);
        }
        updateNeighShifts();
    }

    template<class ScalarType, class PosScalarType>
    CellList<ScalarType, PosScalarType>::CellList(const MDCellType& mdCell, ScalarType cutoff_)
            : lattice(mdCell.getLattice())
            , directPos(mdCell.getPos())
            , cutoff(cutoff_)
            , atomCellMap(mdCell.getNumParticle())
            , neighShifts(3 * 3 * 3) {
        const auto dim = makeGridDim(lattice, cutoff);
        cellGrid.resize(dim);

        mdCell.toDirect(directPos);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Index3D index = posToIndex(i);
            atomCellMap[i] = index;
            cellGrid(index).push_front(i);
        }
        updateNeighShifts();
    }

    template<class ScalarType, class PosScalarType>
    CellList<ScalarType, PosScalarType>& CellList<ScalarType, PosScalarType>::operator=(CellList list) noexcept {
        swap(list);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    void CellList<ScalarType, PosScalarType>::update(const MDCellType& mdCell) {
        directPos = mdCell.getPos();
        mdCell.toDirect(directPos);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Index3D indexNewCell = posToIndex(i);
            if (atomCellMap[i] != indexNewCell) {
                const Index3D indexOldCell = atomCellMap[i];
                cellGrid(indexOldCell).remove(i);
                cellGrid(indexNewCell).push_front(i);
                atomCellMap[i] = indexNewCell;
            }
        }
    }

    template<class ScalarType, class PosScalarType>
    void CellList<ScalarType, PosScalarType>::swap(CellList& list) noexcept {
        cellGrid.swap(list.cellGrid);
        cutoff.swap(list.cutoff);
        atomCellMap.swap(list.atomCellMap);
        directPos.swap(list.directPos);
        neighShifts.swap(list.neighShifts);
    }

    template<class ScalarType, class PosScalarType>
    template<class Functor>
    void CellList<ScalarType, PosScalarType>::forCellInList(Functor func) const {
        for (size_t x = 0; x < cellGrid.getDimX(); ++x)
            for (size_t y = 0; y < cellGrid.getDimY(); ++y)
                for (size_t z = 0; z < cellGrid.getDimZ(); ++z)
                    func(Index3D{x, y, z});

    }

    template<class ScalarType, class PosScalarType>
    template<class Functor>
    void CellList<ScalarType, PosScalarType>::forNeighInRange(Index3D centerCell, Functor func) const {
        auto a1 = lattice.row(0);
        auto a2 = lattice.row(1);
        auto a3 = lattice.row(2);

        Index3D index{};
        const size_t centerX = centerCell[0];
        for (int deltaX = -1; deltaX <= 1; ++deltaX) {
            const int indexShiftX = findNeighbor<0>(centerX, deltaX, index);
            const size_t centerY = centerCell[1];

            for (int deltaY = -1; deltaY <= 1; ++deltaY) {
                const int indexShiftY = findNeighbor<1>(centerY, deltaY, index);
                const size_t centerZ = centerCell[2];

                for (int deltaZ = -1; deltaZ <= 1; ++deltaZ) {
                    if (deltaX == 0 && deltaY == 0 && deltaZ == 0) [[unlikely]]
                        continue;
                    const int indexShiftZ = findNeighbor<2>(centerZ, deltaZ, index);
                    const int indexShift = indexShiftX * 3 * 3 + indexShiftY * 3 + indexShiftZ;
                    func(neighShifts[indexShift], index);
                }
            }
        }
    }

    template<class ScalarType, class PosScalarType>
    template<class Functor>
    void CellList<ScalarType, PosScalarType>::forReducedNeighInRange(Index3D centerCell, Functor func) const {
        auto a1 = lattice.row(0);
        auto a2 = lattice.row(1);
        auto a3 = lattice.row(2);

        Index3D index{};
        const size_t centerX = centerCell[0];
        for (int deltaX = 0; deltaX <= 1; ++deltaX) {
            const int indexShiftX = findNeighbor<0>(centerX, deltaX, index);
            const size_t centerY = centerCell[1];

            for (int deltaY = (deltaX == 0 ? 0 : -1); deltaY <= 1; ++deltaY) {
                const int indexShiftY = findNeighbor<1>(centerY, deltaY, index);
                const size_t centerZ = centerCell[2];

                for (int deltaZ = ((deltaX == 0 && deltaY == 0) ? 1 : -1); deltaZ <= 1; ++deltaZ) {
                    const int indexShiftZ = findNeighbor<2>(centerZ, deltaZ, index);
                    const int indexShift = indexShiftX * 3 * 3 + indexShiftY * 3 + indexShiftZ;
                    func(neighShifts[indexShift], index);
                }
            }
        }
    }

    template<class ScalarType, class PosScalarType>
    typename CellList<ScalarType, PosScalarType>::Index3D
    CellList<ScalarType, PosScalarType>::posToIndex(size_t atomId) const {
        assert(PosScalarType(0) <= directPos(atomId, 0) && directPos(atomId, 0) <= PosScalarType(1));
        assert(PosScalarType(0) <= directPos(atomId, 1) && directPos(atomId, 1) <= PosScalarType(1));
        assert(PosScalarType(0) <= directPos(atomId, 2) && directPos(atomId, 2) <= PosScalarType(1));
        const ScalarType x0 = abs(directPos(atomId, 0) - std::numeric_limits<ScalarType>::epsilon());
        const ScalarType y0 = abs(directPos(atomId, 1) - std::numeric_limits<ScalarType>::epsilon());
        const ScalarType z0 = abs(directPos(atomId, 2) - std::numeric_limits<ScalarType>::epsilon());
        const size_t x = size_t(double(x0 * ScalarType(cellGrid.getDimX())));
        const size_t y = size_t(double(y0 * ScalarType(cellGrid.getDimY())));
        const size_t z = size_t(double(z0 * ScalarType(cellGrid.getDimZ())));
        return {x, y, z};
    }

    template<class ScalarType, class PosScalarType>
    void CellList<ScalarType, PosScalarType>::updateNeighShifts() {
        for (int x = 0; x < 3; ++x) {
            PosScalarType deltaX(x - 1);
            for (int y = 0; y < 3; ++y) {
                PosScalarType deltaY(y - 1);
                for (int z = 0; z < 3; ++z) {
                    PosScalarType deltaZ(z - 1);
                    const int index = x * 3 * 3 + y * 3 + z;
                    neighShifts[index] = lattice.transpose() * Vector3D{deltaX, deltaY, deltaZ};
                }
            }
        }
    }

    template<class ScalarType, class PosScalarType>
    template<size_t DimID>
    int CellList<ScalarType, PosScalarType>::findNeighbor(
            size_t centerIndex,
            int deltaIndex,
            Index3D& neighborIndex) const {
        const size_t dimSize = cellGrid.getDim()[DimID];
        const size_t neigh = static_cast<ssize_t>(centerIndex) + deltaIndex;
        const size_t arr[3]{dimSize - 1, neigh, 0};
        const bool isGood = neigh < dimSize;
        const bool isOverflow = neigh == dimSize;
        const int index = isGood + isOverflow * 2;
        neighborIndex[DimID] = arr[index];
        return index;
    }

    template<class ScalarType, class PosScalarType>
    typename CellList<ScalarType, PosScalarType>::Index3D
    CellList<ScalarType, PosScalarType>::makeGridDim(const LatticeMatrix& lattice, ScalarType cutoff) {
        size_t dimX, dimY, dimZ;
        const ReciprocalCell repCell(lattice);
        const auto& repLatt = repCell.getLattice();
        const ScalarType factor = reciprocal(cutoff) * (2 * M_PI);
        dimX = static_cast<size_t>(double(factor * reciprocal(repLatt.row(0).norm())));
        dimY = static_cast<size_t>(double(factor * reciprocal(repLatt.row(1).norm())));
        dimZ = static_cast<size_t>(double(factor * reciprocal(repLatt.row(2).norm())));
        if (dimX == 0 || dimY == 0 || dimZ == 0)
            throw std::invalid_argument("[Error]: Cell is too small that self interaction is calculated");
        return {dimX, dimY, dimZ};
    }
}
