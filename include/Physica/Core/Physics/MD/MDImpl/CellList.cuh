/*
 * Copyright 2023 WeiBo He.
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

#include "CellList.h"
#include "Physica/Utils/CUDA/device_obj.cuh"

namespace Physica::Core {
    template<class ScalarType>
    class device_obj<CellList<ScalarType>> {
        using host_obj = CellList<ScalarType>;
        using This = device_obj<host_obj>;
        using Index3D = typename GridBase::Index3D;
        using Vector3D = typename host_obj::Vector3D;
        using LatticeMatrix = device_obj<typename host_obj::LatticeMatrix>;
        using DeviceIndexArray = Utils::device_obj<Utils::Array<size_t>>;
        using DeviceNeighShift = Utils::device_obj<Utils::Array<Vector3D>>;

        LatticeMatrix lattice;
        Index3D cellGridDim;
        DeviceIndexArray cellAtomMap;
        DeviceIndexArray cellStartOffset;
        DeviceNeighShift neighShifts;
    public:
        device_obj() = default;
        device_obj(const host_obj& obj);
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        device_obj& operator=(device_obj obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class Functor>
        __device__ void forCellInList(Functor func) const;
        template<class Functor>
        __device__ void forNeighInRange(Index3D centerCell, Functor func) const;
        template<class Functor>
        __device__ inline void forAtomInCell(Index3D cellIndex, Functor func) const;
        void swap(device_obj& obj) noexcept;
        /* Getters */
        [[nodiscard]] __device__ const LatticeMatrix& getLattice() const noexcept { return lattice; }
        [[nodiscard]] __host__ __device__ Index3D getCellGridDim() const noexcept { return cellGridDim; }
        [[nodiscard]] __host__ __device__ size_t getCellGridDimX() const noexcept { return cellGridDim[0]; }
        [[nodiscard]] __host__ __device__ size_t getCellGridDimY() const noexcept { return cellGridDim[1]; }
        [[nodiscard]] __host__ __device__ size_t getCellGridDimZ() const noexcept { return cellGridDim[2]; }
        [[nodiscard]] __device__ const DeviceIndexArray& getCellAtomMap() const noexcept { return cellAtomMap; }
        [[nodiscard]] __device__ const DeviceIndexArray& getCellStartOffset() const noexcept { return cellStartOffset; }
        [[nodiscard]] __device__ size_t getNumAtomInCell(size_t cell) { return cellStartOffset[cell + 1] - cellStartOffset[cell]; }
        [[nodiscard]] __host__ __device__ size_t getNumCell() const noexcept { return cellGridDim[0] * cellGridDim[1] * cellGridDim[2]; }
        [[nodiscard]] __device__ const DeviceNeighShift& getNeighShifts() const noexcept { return neighShifts; }
        /* Setters */
        void setCellGridDim(Index3D cellGridDim_) { cellGridDim = std::move(cellGridDim_); }
    private:
        friend class CellList<ScalarType>;
    };

    template<class ScalarType>
    device_obj<CellList<ScalarType>>::device_obj(const host_obj& obj) {
        obj.toDevice(*this);
    }

    template<class ScalarType>
    template<class Functor>
    __device__ void device_obj<CellList<ScalarType>>::forCellInList(Functor func) const {
        for (size_t x = 0; x < getCellGridDimX(); ++x)
            for (size_t y = 0; y < getCellGridDimY(); ++y)
                for (size_t z = 0; z < getCellGridDimZ(); ++z)
                    func(Index3D{x, y, z});
    }

    template<class ScalarType>
    template<class Functor>
    __device__ void device_obj<CellList<ScalarType>>::forNeighInRange(Index3D centerCell, Functor func) const {
        auto a1 = lattice.row(0);
        auto a2 = lattice.row(1);
        auto a3 = lattice.row(2);

        Index3D index{};
        const size_t centerX = centerCell[0];
        for (int deltaX = -1; deltaX <= 1; ++deltaX) {
            const int indexShiftX = host_obj::template findNeighbor<0>(cellGridDim, centerX, deltaX, index);
            const size_t centerY = centerCell[1];

            for (int deltaY = -1; deltaY <= 1; ++deltaY) {
                const int indexShiftY = host_obj::template findNeighbor<1>(cellGridDim, centerY, deltaY, index);
                const size_t centerZ = centerCell[2];

                for (int deltaZ = -1; deltaZ <= 1; ++deltaZ) {
                    if (deltaX == 0 && deltaY == 0 && deltaZ == 0) [[unlikely]]
                        continue;
                    const int indexShiftZ = host_obj::template findNeighbor<2>(cellGridDim, centerZ, deltaZ, index);
                    const int indexShift = indexShiftX * 3 * 3 + indexShiftY * 3 + indexShiftZ;
                    func(neighShifts[indexShift], index);
                }
            }
        }
    }

    template<class ScalarType>
    template<class Functor>
    __device__ inline void device_obj<CellList<ScalarType>>::forAtomInCell(Index3D cellIndex, Functor func) const {
        const auto index1D = PeriodIndex3D(cellIndex, cellGridDim).toIndex1D();
        const size_t cellBegin = cellStartOffset[index1D];
        const size_t cellEnd = cellStartOffset[index1D + 1];
        for (size_t i = cellBegin; i < cellEnd; ++i)
            func(cellAtomMap[i]);
    }

    template<class ScalarType>
    void device_obj<CellList<ScalarType>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        lattice.swap(obj.lattice);
        cellGridDim.swap(obj.cellGridDim);
        cellAtomMap.swap(obj.cellAtomMap);
        cellStartOffset.swap(obj.cellStartOffset);
        neighShifts.swap(obj.neighShifts);
    }

    template<class ScalarType>
    inline device_obj<CellList<ScalarType>> CellList<ScalarType>::toDevice() const {
        return device_obj<CellList<ScalarType>>(*this);
    }

    template<class ScalarType>
    void CellList<ScalarType>::toDevice(device_obj<CellList<ScalarType>>& obj) const {
        lattice.toDevice(obj.lattice);
        obj.cellGridDim = cellGridDim;
        cellAtomMap.toDevice(obj.cellAtomMap);
        cellStartOffset.toDevice(obj.cellStartOffset);
        neighShifts.toDevice(obj.neighShifts);
    }
}
