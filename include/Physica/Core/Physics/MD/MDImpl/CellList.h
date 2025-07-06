/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Tensor/PeriodIndex3D.h"
#include "Physica/Core/Physics/MD/MDCell.h"
#include "Physica/Core/Utils/Container/ArrayND.h"

namespace Physica {
    template<Scalar T>
    class CellList {
        static_assert(!Diffable<T>, "[Error]: CellList is not diffable");
        using This = CellList<T>;
    public:
        using MDCellType = MDCell<T>;
        using LatticeMatrix = MDCellType::LatticeMatrix;
        using PositionMatrix = MDCellType::PositionMatrix;
        using CellGrid = ArrayND<std::forward_list<size_t>, 3>;
    private:
        LatticeMatrix lattice;
        PositionMatrix directPos;
        T cutoff;
        Index3D cellGridDim;
        Array<size_t> cellAtomMap;
        Array<size_t> cellStartOffset;
        Array<Index3D> atomCellMap;
        Array<Vector3D<T>> neighShifts;
    public:
        CellList(LatticeMatrix lattice_, PositionMatrix cartesianPos, T cutoff_);
        CellList(const MDCellType& mdCell, T cutoff_);
        CellList(const CellList&) = default;
        CellList(CellList&&) noexcept = default;
        ~CellList() = default;
        /* Operators */
        CellList& operator=(CellList obj) noexcept;
        /* Operations */
        void update(PositionMatrix pos);
        void update(const MDCellType& mdCell);
        void forCellInList(std::invocable<Index3D> auto fn) const;
        void forNeighInRange(Index3D centerCell, std::invocable<const Vector3D<T>&, Index3D> auto fn) const;
        void forReducedNeighInRange(Index3D centerCell, std::invocable<const Vector3D<T>&, Index3D> auto fn) const;
        void forAtomInCell(Index3D cellIndex, std::invocable<size_t> auto fn) const;

        [[nodiscard]] size_t calcMaxNumAtomInCell() const noexcept;

        [[nodiscard]] inline device_obj<This> toDevice() const;
        void toDevice(device_obj<This>& obj) const;
        void swap(CellList& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] size_t getNumParticle() const noexcept { return directPos.getRow(); }
        [[nodiscard]] T getCutoff() const noexcept { return cutoff; }
        [[nodiscard]] Index3D getCellGridDim() const noexcept { return cellGridDim; }
        [[nodiscard]] size_t getCellGridDimX() const noexcept { return cellGridDim[0]; }
        [[nodiscard]] size_t getCellGridDimY() const noexcept { return cellGridDim[1]; }
        [[nodiscard]] size_t getCellGridDimZ() const noexcept { return cellGridDim[2]; }
        [[nodiscard]] const Array<size_t>& getCellAtomMap() const noexcept { return cellAtomMap; }
        [[nodiscard]] const Array<size_t>& getCellStartOffset() const noexcept { return cellStartOffset; }
        [[nodiscard]] inline size_t getNumAtomInCell(size_t cell) const;
        [[nodiscard]] size_t getNumCell() const noexcept { return cellStartOffset.getLength() - 1; }
        [[nodiscard]] const Array<Index3D>& getAtomCellMap() const noexcept { return atomCellMap; }
    private:
        Index3D posToIndex(size_t atomId) const;
        /* Static members */
        [[nodiscard]] static Array<Vector3D<T>> makeNeighShifts(const LatticeMatrix& lattice);
        template<size_t DimID>
        [[nodiscard]] __host__ __device__ inline static int findNeighbor(
                const Index3D& cellGridDim, size_t centerIndex, int deltaIndex, Index3D& neighborIndex);
        [[nodiscard]] static Index3D makeGridDim(const LatticeMatrix& lattice, T cutoff);
        /* Friends */
        friend class device_obj<This>;
    };

    template<Scalar T>
    CellList<T>::CellList(LatticeMatrix lattice_, PositionMatrix cartesianPos, T cutoff_)
            : lattice(std::move(lattice_))
            , directPos(std::move(cartesianPos))
            , cutoff(cutoff_) {
        assert(cutoff.isPositive() && "[Error]: Invalid cutoff");
        cellGridDim = makeGridDim(lattice, cutoff);
        cellAtomMap.resize(getNumParticle());
        cellStartOffset.resize(cellGridDim[0] * cellGridDim[1] * cellGridDim[2] + 1);
        (*cellStartOffset.rbegin()) = getNumParticle();
        atomCellMap.resize(getNumParticle());

        PeriodicCell<T, 3>::toDirect(directPos, lattice);
        auto cellGrid = CellGrid(cellGridDim);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Index3D index = posToIndex(i);
            atomCellMap[i] = index;
            cellGrid(index).push_front(i);
        }
        {
            size_t cellAtomMapIndex = 0;
            auto& cells = cellGrid.asArray();
            for (size_t i = 0; i < cells.getLength(); ++i) {
                cellStartOffset[i] = cellAtomMapIndex;
                const auto& list = cells[i];
                for (auto atom : list) {
                    cellAtomMap[cellAtomMapIndex] = atom;
                    cellAtomMapIndex += 1;
                }
            }
        }
        neighShifts = makeNeighShifts(lattice);
    }

    template<Scalar T>
    CellList<T>::CellList(const MDCellType& mdCell, T cutoff_)
            : CellList(mdCell.getLattice(), mdCell.getPos(), std::move(cutoff_)) {}

    template<Scalar T>
    CellList<T>& CellList<T>::operator=(CellList obj) noexcept {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void CellList<T>::update(PositionMatrix pos) {
        directPos = std::move(pos);
        PeriodicCell<T, 3>::toDirect(directPos, lattice);
        CellGrid cellGrid(cellGridDim);
        for (size_t i = 0; i < getNumParticle(); ++i) {
            const Index3D index = posToIndex(i);
            atomCellMap[i] = index;
            cellGrid(index).push_front(i);
        }
        {
            size_t cellAtomMapIndex = 0;
            auto& cells = cellGrid.asArray();
            for (size_t i = 0; i < cells.getLength(); ++i) {
                cellStartOffset[i] = cellAtomMapIndex;
                const auto& list = cells[i];
                for (auto atom : list) {
                    cellAtomMap[cellAtomMapIndex] = atom;
                    cellAtomMapIndex += 1;
                }
            }
        }
    }

    template<Scalar T>
    inline void CellList<T>::update(const MDCellType& mdCell) {
        update(mdCell.getPos());
    }

    template<Scalar T>
    void CellList<T>::forCellInList(std::invocable<Index3D> auto fn) const {
        for (size_t x = 0; x < getCellGridDimX(); ++x)
            for (size_t y = 0; y < getCellGridDimY(); ++y)
                for (size_t z = 0; z < getCellGridDimZ(); ++z)
                    fn(Index3D{x, y, z});
    }

    template<Scalar T>
    void CellList<T>::forNeighInRange(Index3D centerCell, std::invocable<const Vector3D<T>&, Index3D> auto fn) const {
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
                    fn(neighShifts[indexShift], index);
                }
            }
        }
    }

    template<Scalar T>
    void CellList<T>::forReducedNeighInRange(Index3D centerCell, std::invocable<const Vector3D<T>&, Index3D> auto fn) const {
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
                    fn(neighShifts[indexShift], index);
                }
            }
        }
    }

    template<Scalar T>
    void CellList<T>::forAtomInCell(Index3D cellIndex, std::invocable<size_t> auto fn) const {
        const auto index1D = PeriodIndex3D(cellIndex, cellGridDim).toIndex1D();
        const size_t cellBegin = cellStartOffset[index1D];
        const size_t cellEnd = cellStartOffset[index1D + 1];
        for (size_t i = cellBegin; i < cellEnd; ++i)
            fn(cellAtomMap[i]);
    }

    template<Scalar T>
    size_t CellList<T>::calcMaxNumAtomInCell() const noexcept {
        size_t result = 0;
        for (size_t cell = 0; cell < getNumCell(); ++cell) {
            const size_t numAtom = getNumAtomInCell(cell);
            result = std::max(result, numAtom);
        }
        return result;
    }

    template<Scalar T>
    void CellList<T>::swap(CellList& __restrict obj) noexcept {
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

    template<Scalar T>
    inline size_t CellList<T>::getNumAtomInCell(size_t cell) const {
        assert(cell < getNumCell() && "[Error]: The cell is not exist");
        return cellStartOffset[cell + 1] - cellStartOffset[cell];
    }

    template<Scalar T>
    Index3D CellList<T>::posToIndex(size_t atomId) const {
        const T x0 = abs(directPos(atomId, 0) - T(std::numeric_limits<T>::epsilon()));
        const T y0 = abs(directPos(atomId, 1) - T(std::numeric_limits<T>::epsilon()));
        const T z0 = abs(directPos(atomId, 2) - T(std::numeric_limits<T>::epsilon()));
        assert(T(0) <= x0 && x0 <= T(1));
        assert(T(0) <= y0 && y0 <= T(1));
        assert(T(0) <= z0 && z0 <= T(1));
        const size_t x = size_t(double(x0 * T(getCellGridDimX())));
        const size_t y = size_t(double(y0 * T(getCellGridDimY())));
        const size_t z = size_t(double(z0 * T(getCellGridDimZ())));
        return {x, y, z};
    }

    template<Scalar T>
    Array<Vector3D<T>> CellList<T>::makeNeighShifts(const LatticeMatrix& lattice) {
        Array<Vector3D<T>> neighShifts(3 * 3 * 3);
        for (int x = 0; x < 3; ++x) {
            T deltaX(x - 1);
            for (int y = 0; y < 3; ++y) {
                T deltaY(y - 1);
                for (int z = 0; z < 3; ++z) {
                    T deltaZ(z - 1);
                    const int index = x * 3 * 3 + y * 3 + z;
                    neighShifts[index] = lattice.transpose() * Vector3D<T>{deltaX, deltaY, deltaZ};
                }
            }
        }
        return neighShifts;
    }

    template<Scalar T>
    template<size_t DimID>
    __host__ __device__ inline int CellList<T>::findNeighbor(
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

    template<Scalar T>
    Index3D CellList<T>::makeGridDim(const LatticeMatrix& lattice, T cutoff) {
        const auto repLatt = MDCellType::makeRepLattice(lattice);
        const T factor = reciprocal(cutoff) * T(2 * M_PI);
        size_t dimX = static_cast<size_t>(double(factor * reciprocal(repLatt.row(0).norm())));
        size_t dimY = static_cast<size_t>(double(factor * reciprocal(repLatt.row(1).norm())));
        size_t dimZ = static_cast<size_t>(double(factor * reciprocal(repLatt.row(2).norm())));
        if (dimX == 0 || dimY == 0 || dimZ == 0)
            throw std::invalid_argument("[Error]: Cell is too small that self interaction is calculated");
        return {dimX, dimY, dimZ};
    }
}
