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
#include "Physica/Core/Physics/Container/RSpaceGrid.h"
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    template<class ScalarType, class PosScalarType>
    class CellList {
        using MDCellType = MDCell<ScalarType, PosScalarType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using CellGrid = RSpaceGrid<std::forward_list<size_t>>;
        using Vector3D = Vector<PosScalarType, 3>;
    public:
        using Index3D = typename CellGrid::Index3D;
    private:
        LatticeMatrix lattice;
        PositionMatrix buffer;
        CellGrid cellGrid;
        ScalarType cutoff;
        Utils::Array<size_t> atomCellMap;
    public:
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
        /* Static members */
        template<class Functor>
        void forCellInList(Functor func) const;
        template<class Functor>
        void forNeighInRange(Index3D centerCell, Functor func) const;
    private:
        Index3D posToIndex(size_t atomId) const;
        /* Helpers */
        inline Vector3D findNeighbor(size_t dimId, Index3D centerIndex, int deltaIndex, Index3D& neighborIndex) const;
        /* Static members */
        static Index3D makeGridDim(const MDCellType& mdCell, ScalarType cutoff);
    };

    template<class ScalarType, class PosScalarType>
    CellList<ScalarType, PosScalarType>::CellList(const MDCellType& mdCell, ScalarType cutoff_)
            : lattice(mdCell.getLattice())
            , buffer(mdCell.getPos())
            , cutoff(cutoff_)
            , atomCellMap(mdCell.getNumParticle()) {
        const auto dim = makeGridDim(mdCell, cutoff);
        cellGrid.resize(dim[0], dim[1], dim[2]);

        mdCell.toDirect(buffer);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Index3D index3d = posToIndex(i);
            const size_t index1d = cellGrid.index3dTo1d(index3d);
            atomCellMap[i] = index1d;
            cellGrid.asVector()[index1d].push_front(i);
        }
    }

    template<class ScalarType, class PosScalarType>
    CellList<ScalarType, PosScalarType>& CellList<ScalarType, PosScalarType>::operator=(CellList list) noexcept {
        swap(list);
        return *this;
    }

    template<class ScalarType, class PosScalarType>
    void CellList<ScalarType, PosScalarType>::update(const MDCellType& mdCell) {
        buffer = mdCell.getPos();
        mdCell.toDirect(buffer);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Index3D index3d = posToIndex(i);
            const size_t indexNewCell = cellGrid.index3dTo1d(index3d);
            if (atomCellMap[i] != indexNewCell) {
                const size_t indexOldCell = atomCellMap[i];
                cellGrid.asVector()[indexOldCell].remove(i);
                cellGrid.asVector()[indexNewCell].push_front(i);
                atomCellMap[i] = indexNewCell;
            }
        }
    }

    template<class ScalarType, class PosScalarType>
    void CellList<ScalarType, PosScalarType>::swap(CellList& list) noexcept {
        cellGrid.swap(list.cellGrid);
        cutoff.swap(list.cutoff);
        atomCellMap.swap(list.atomCellMap);
        buffer.swap(list.buffer);
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

        Vector3D v1, v2, v3;
        Index3D index{};
        for (int deltaX = -1; deltaX <= 1; ++deltaX) {
            v1 = findNeighbor(0, centerCell, deltaX, index);
            for (int deltaY = -1; deltaY <= 1; ++deltaY) {
                v2 = v1 + findNeighbor(1, centerCell, deltaY, index);
                for (int deltaZ = -1; deltaZ <= 1; ++deltaZ) {
                    v3 = v2 + findNeighbor(2, centerCell, deltaZ, index);
                    func(v3, index);
                }
            }
        }
    }

    template<class ScalarType, class PosScalarType>
    typename CellList<ScalarType, PosScalarType>::Index3D
    CellList<ScalarType, PosScalarType>::posToIndex(size_t atomId) const {
        assert(ScalarType(0) <= buffer(atomId, 0) && buffer(atomId, 0) <= ScalarType(1));
        assert(ScalarType(0) <= buffer(atomId, 1) && buffer(atomId, 1) <= ScalarType(1));
        assert(ScalarType(0) <= buffer(atomId, 2) && buffer(atomId, 2) <= ScalarType(1));
        const size_t x = size_t(double(buffer(atomId, 0) * ScalarType(cellGrid.getDimX())));
        const size_t y = size_t(double(buffer(atomId, 1) * ScalarType(cellGrid.getDimY())));
        const size_t z = size_t(double(buffer(atomId, 2) * ScalarType(cellGrid.getDimZ())));
        return {x, y, z};
    }

    template<class ScalarType, class PosScalarType>
    typename CellList<ScalarType, PosScalarType>::Vector3D
    CellList<ScalarType, PosScalarType>::findNeighbor(
            size_t dimId,
            Index3D centerIndex,
            int deltaIndex,
            Index3D& neighborIndex) const {
        const size_t dimSize = cellGrid.getDim()[dimId];
        const auto lattVec = lattice.row(dimId);

        const size_t neigh = static_cast<ssize_t>(centerIndex[dimId]) + deltaIndex;
        Vector3D result{};

        const bool isGood = neigh < dimSize;
        if (isGood) {
            result = PosScalarType(0);
            neighborIndex[dimId] = neigh;
            return result;
        }
        const bool isOverflow = neigh == dimSize;
        if (isOverflow) {
            result = lattVec;
            neighborIndex[dimId] = 0;
            return result;
        }
        /* Underflow */ {
            result = -lattVec.asVector();
            neighborIndex[dimId] = dimSize - 1;
            return result;
        }
    }

    template<class ScalarType, class PosScalarType>
    typename CellList<ScalarType, PosScalarType>::Index3D
    CellList<ScalarType, PosScalarType>::makeGridDim(const MDCellType& mdCell, ScalarType cutoff) {
        size_t dimX, dimY, dimZ;
        const ReciprocalCell repCell(mdCell.getLattice());
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
