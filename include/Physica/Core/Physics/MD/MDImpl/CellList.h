/*
 * Copyright 2022-2023 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/PeriodIndex3D.h"
#include "Physica/Core/Physics/MD/MDCell.h"

namespace Physica::Core {
    template<class ScalarType>
    class CellList {
        using This = CellList<ScalarType>;
        using PlainScalar = typename ScalarType::PlainScalar;
    public:
        using MDCellType = MDCell<ScalarType>;
        using LatticeMatrix = typename MDCellType::LatticeMatrix;
        using PositionMatrix = typename MDCellType::PositionMatrix;
        using CellGrid = GridStorage<std::forward_list<size_t>>;
        using Vector3D = Vector<ScalarType, 3>;
        using Index3D = typename GridBase::Index3D;
    private:
        LatticeMatrix lattice;
        PositionMatrix directPos;
        PlainScalar cutoff;
        Index3D cellGridDim;
        Utils::Array<size_t> cellAtomMap;
        Utils::Array<size_t> cellStartOffset;
        Utils::Array<Index3D> atomCellMap;
        Utils::Array<Vector3D> neighShifts;
    public:
        CellList(LatticeMatrix lattice_, PositionMatrix cartesianPos, PlainScalar cutoff_);
        CellList(const MDCellType& mdCell, PlainScalar cutoff_);
        CellList(const CellList&) = default;
        CellList(CellList&&) noexcept = default;
        ~CellList() = default;
        /* Operators */
        CellList& operator=(CellList obj) noexcept;
        /* Operations */
        void update(PositionMatrix pos);
        inline void update(const MDCellType& mdCell);
        template<class Functor>
        void forCellInList(Functor func) const;
        template<class Functor>
        void forNeighInRange(Index3D centerCell, Functor func) const;
        template<class Functor>
        void forReducedNeighInRange(Index3D centerCell, Functor func) const;
        template<class Functor>
        inline void forAtomInCell(Index3D cellIndex, Functor func) const;

        [[nodiscard]] size_t calcMaxNumAtomInCell() const noexcept;

        [[nodiscard]] inline device_obj<This> toDevice() const;
        void toDevice(device_obj<This>& obj) const;
        void swap(CellList& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return directPos.getRow(); }
        [[nodiscard]] PlainScalar getCutoff() const noexcept { return cutoff; }
        [[nodiscard]] Index3D getCellGridDim() const noexcept { return cellGridDim; }
        [[nodiscard]] size_t getCellGridDimX() const noexcept { return cellGridDim[0]; }
        [[nodiscard]] size_t getCellGridDimY() const noexcept { return cellGridDim[1]; }
        [[nodiscard]] size_t getCellGridDimZ() const noexcept { return cellGridDim[2]; }
        [[nodiscard]] const Utils::Array<size_t>& getCellAtomMap() const noexcept { return cellAtomMap; }
        [[nodiscard]] const Utils::Array<size_t>& getCellStartOffset() const noexcept { return cellStartOffset; }
        [[nodiscard]] inline size_t getNumAtomInCell(size_t cell) const;
        [[nodiscard]] size_t getNumCell() const noexcept { return cellStartOffset.getLength() - 1; }
        [[nodiscard]] const Utils::Array<Index3D>& getAtomCellMap() const noexcept { return atomCellMap; }
    private:
        Index3D posToIndex(size_t atomId) const;
        /* Static members */
        [[nodiscard]] static Utils::Array<Vector3D> makeNeighShifts(const LatticeMatrix& lattice);
        template<size_t DimID>
        [[nodiscard]] __host__ __device__ inline static int findNeighbor(
                const Index3D& cellGridDim, size_t centerIndex, int deltaIndex, Index3D& neighborIndex);
        [[nodiscard]] static Index3D makeGridDim(const LatticeMatrix& lattice, PlainScalar cutoff);
        /* Friends */
        friend class device_obj<This>;
    };

    template<class ScalarType>
    CellList<ScalarType>::CellList(LatticeMatrix lattice_, PositionMatrix cartesianPos, PlainScalar cutoff_)
            : lattice(std::move(lattice_))
            , directPos(std::move(cartesianPos))
            , cutoff(cutoff_) {
        assert(cutoff.isPositive() && "[Error]: Invalid cutoff");
        cellGridDim = makeGridDim(lattice, cutoff);
        cellAtomMap.resize(getNumParticle());
        cellStartOffset.resize(cellGridDim[0] * cellGridDim[1] * cellGridDim[2] + 1);
        (*cellStartOffset.rbegin()) = getNumParticle();
        atomCellMap.resize(getNumParticle());

        PeriodicCell<ScalarType, 3>::toDirect(directPos, lattice);
        CellGrid cellGrid(cellGridDim);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Index3D index = posToIndex(i);
            atomCellMap[i] = index;
            cellGrid(index).push_front(i);
        }
        {
            size_t cellAtomMapIndex = 0;
            cellGrid.forIndexInGrid([this, &cellGrid, &cellAtomMapIndex](Index3D index) {
                const auto index1D = PeriodIndex3D(index, cellGridDim).toIndex1D();
                cellStartOffset[index1D] = cellAtomMapIndex;
                const auto& list = cellGrid(index);
                for (auto atom : list) {
                    cellAtomMap[cellAtomMapIndex] = atom;
                    cellAtomMapIndex += 1;
                }
            });
        }
        neighShifts = makeNeighShifts(lattice);
    }

    template<class ScalarType>
    CellList<ScalarType>::CellList(const MDCellType& mdCell, PlainScalar cutoff_)
            : CellList(mdCell.getLattice(), mdCell.getPos(), std::move(cutoff_)) {}

    template<class ScalarType>
    CellList<ScalarType>& CellList<ScalarType>::operator=(CellList obj) noexcept {
        swap(obj);
        return *this;
    }

    template<class ScalarType>
    void CellList<ScalarType>::update(PositionMatrix pos) {
        directPos = std::move(pos);
        PeriodicCell<ScalarType, 3>::toDirect(directPos, lattice);
        CellGrid cellGrid(cellGridDim);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Index3D index = posToIndex(i);
            atomCellMap[i] = index;
            cellGrid(index).push_front(i);
        }
        {
            size_t cellAtomMapIndex = 0;
            cellGrid.forIndexInGrid([this, &cellGrid, &cellAtomMapIndex](Index3D index) {
                const auto index1D = PeriodIndex3D(index, cellGridDim).toIndex1D();
                cellStartOffset[index1D] = cellAtomMapIndex;
                const auto& list = cellGrid(index);
                for (auto atom : list) {
                    cellAtomMap[cellAtomMapIndex] = atom;
                    cellAtomMapIndex += 1;
                }
            });
        }
    }

    template<class ScalarType>
    inline void CellList<ScalarType>::update(const MDCellType& mdCell) {
        update(mdCell.getPos());
    }

    template<class ScalarType>
    template<class Functor>
    void CellList<ScalarType>::forCellInList(Functor func) const {
        for (size_t x = 0; x < getCellGridDimX(); ++x)
            for (size_t y = 0; y < getCellGridDimY(); ++y)
                for (size_t z = 0; z < getCellGridDimZ(); ++z)
                    func(Index3D{x, y, z});
    }

    template<class ScalarType>
    template<class Functor>
    void CellList<ScalarType>::forNeighInRange(Index3D centerCell, Functor func) const {
        auto a1 = lattice.row(0);
        auto a2 = lattice.row(1);
        auto a3 = lattice.row(2);

        Index3D index{};
        const size_t centerX = centerCell[0];
        for (int deltaX = -1; deltaX <= 1; ++deltaX) {
            const int indexShiftX = findNeighbor<0>(cellGridDim, centerX, deltaX, index);
            const size_t centerY = centerCell[1];

            for (int deltaY = -1; deltaY <= 1; ++deltaY) {
                const int indexShiftY = findNeighbor<1>(cellGridDim, centerY, deltaY, index);
                const size_t centerZ = centerCell[2];

                for (int deltaZ = -1; deltaZ <= 1; ++deltaZ) {
                    if (deltaX == 0 && deltaY == 0 && deltaZ == 0) [[unlikely]]
                        continue;
                    const int indexShiftZ = findNeighbor<2>(cellGridDim, centerZ, deltaZ, index);
                    const int indexShift = indexShiftX * 3 * 3 + indexShiftY * 3 + indexShiftZ;
                    func(neighShifts[indexShift], index);
                }
            }
        }
    }

    template<class ScalarType>
    template<class Functor>
    void CellList<ScalarType>::forReducedNeighInRange(Index3D centerCell, Functor func) const {
        auto a1 = lattice.row(0);
        auto a2 = lattice.row(1);
        auto a3 = lattice.row(2);

        Index3D index{};
        const size_t centerX = centerCell[0];
        for (int deltaX = 0; deltaX <= 1; ++deltaX) {
            const int indexShiftX = findNeighbor<0>(cellGridDim, centerX, deltaX, index);
            const size_t centerY = centerCell[1];

            for (int deltaY = (deltaX == 0 ? 0 : -1); deltaY <= 1; ++deltaY) {
                const int indexShiftY = findNeighbor<1>(cellGridDim, centerY, deltaY, index);
                const size_t centerZ = centerCell[2];

                for (int deltaZ = ((deltaX == 0 && deltaY == 0) ? 1 : -1); deltaZ <= 1; ++deltaZ) {
                    const int indexShiftZ = findNeighbor<2>(cellGridDim, centerZ, deltaZ, index);
                    const int indexShift = indexShiftX * 3 * 3 + indexShiftY * 3 + indexShiftZ;
                    func(neighShifts[indexShift], index);
                }
            }
        }
    }

    template<class ScalarType>
    template<class Functor>
    inline void CellList<ScalarType>::forAtomInCell(Index3D cellIndex, Functor func) const {
        const auto index1D = PeriodIndex3D(cellIndex, cellGridDim).toIndex1D();
        const size_t cellBegin = cellStartOffset[index1D];
        const size_t cellEnd = cellStartOffset[index1D + 1];
        for (size_t i = cellBegin; i < cellEnd; ++i)
            func(cellAtomMap[i]);
    }

    template<class ScalarType>
    size_t CellList<ScalarType>::calcMaxNumAtomInCell() const noexcept {
        size_t result = 0;
        for (size_t cell = 0; cell < getNumCell(); ++cell) {
            const size_t numAtom = getNumAtomInCell(cell);
            result = std::max(result, numAtom);
        }
        return result;
    }

    template<class ScalarType>
    void CellList<ScalarType>::swap(CellList& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lattice.swap(obj.lattice);
        directPos.swap(obj.directPos);
        cutoff.swap(obj.cutoff);
        cellGridDim.swap(obj.cellGridDim);
        cellAtomMap.swap(obj.cellAtomMap);
        cellStartOffset.swap(obj.cellStartOffset);
        atomCellMap.swap(obj.atomCellMap);
        neighShifts.swap(obj.neighShifts);
    }

    template<class ScalarType>
    inline size_t CellList<ScalarType>::getNumAtomInCell(size_t cell) const {
        assert(cell < getNumCell() && "[Error]: The cell is not exist");
        return cellStartOffset[cell + 1] - cellStartOffset[cell];
    }

    template<class ScalarType>
    typename CellList<ScalarType>::Index3D
    CellList<ScalarType>::posToIndex(size_t atomId) const {
        const ScalarType x0 = abs(directPos(atomId, 0) - PlainScalar(std::numeric_limits<ScalarType>::epsilon()));
        const ScalarType y0 = abs(directPos(atomId, 1) - PlainScalar(std::numeric_limits<ScalarType>::epsilon()));
        const ScalarType z0 = abs(directPos(atomId, 2) - PlainScalar(std::numeric_limits<ScalarType>::epsilon()));
        assert(ScalarType(0) <= x0 && x0 <= ScalarType(1));
        assert(ScalarType(0) <= y0 && y0 <= ScalarType(1));
        assert(ScalarType(0) <= z0 && z0 <= ScalarType(1));
        const size_t x = size_t(double(x0 * ScalarType(getCellGridDimX())));
        const size_t y = size_t(double(y0 * ScalarType(getCellGridDimY())));
        const size_t z = size_t(double(z0 * ScalarType(getCellGridDimZ())));
        return {x, y, z};
    }

    template<class ScalarType>
    Utils::Array<typename CellList<ScalarType>::Vector3D>
    CellList<ScalarType>::makeNeighShifts(const LatticeMatrix& lattice) {
        Utils::Array<Vector3D> neighShifts(3 * 3 * 3);
        for (int x = 0; x < 3; ++x) {
            ScalarType deltaX(x - 1);
            for (int y = 0; y < 3; ++y) {
                ScalarType deltaY(y - 1);
                for (int z = 0; z < 3; ++z) {
                    ScalarType deltaZ(z - 1);
                    const int index = x * 3 * 3 + y * 3 + z;
                    neighShifts[index] = lattice.transpose() * Vector3D{deltaX, deltaY, deltaZ};
                }
            }
        }
        return neighShifts;
    }

    template<class ScalarType>
    template<size_t DimID>
    __host__ __device__ inline int CellList<ScalarType>::findNeighbor(
            const Index3D& cellGridDim,
            size_t centerIndex,
            int deltaIndex,
            Index3D& neighborIndex) {
        const size_t dimSize = cellGridDim[DimID];
        const size_t neigh = static_cast<ssize_t>(centerIndex) + deltaIndex;
        const size_t arr[3]{dimSize - 1, neigh, 0};
        const bool isGood = neigh < dimSize;
        const bool isOverflow = neigh == dimSize;
        const int index = isGood + isOverflow * 2;
        neighborIndex[DimID] = arr[index];
        return index;
    }

    template<class ScalarType>
    typename CellList<ScalarType>::Index3D
    CellList<ScalarType>::makeGridDim(const LatticeMatrix& lattice, PlainScalar cutoff) {
        const auto repLatt = MDCellType::makeRepLattice(lattice);
        const PlainScalar factor = reciprocal(cutoff) * PlainScalar(2 * M_PI);
        size_t dimX = static_cast<size_t>(double(factor * reciprocal(repLatt.row(0).norm())));
        size_t dimY = static_cast<size_t>(double(factor * reciprocal(repLatt.row(1).norm())));
        size_t dimZ = static_cast<size_t>(double(factor * reciprocal(repLatt.row(2).norm())));
        if (dimX == 0 || dimY == 0 || dimZ == 0)
            throw std::invalid_argument("[Error]: Cell is too small that self interaction is calculated");
        return {dimX, dimY, dimZ};
    }
}
