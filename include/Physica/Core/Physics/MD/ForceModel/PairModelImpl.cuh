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

#include <cub/cub.cuh>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/PeriodIndex3D.h"
#include "Physica/Core/Parallel/Executor/CudaExecutor.cuh"
#include "Physica/Utils/CUDA/DeviceProp.cuh"
#include "Physica/Utils/CUDA/PlainStruct.h"

namespace Physica::Core {
    namespace Internal {
        template<class Derived, bool IsSmallCell>
        __global__ void PairModel_forceKernel(Physica::PlainStruct<Core::device_obj<PairModel<Derived>>> pair) {
            pair.getDerived().template forceKernelImpl<IsSmallCell>();
        }

        template<class Derived>
        __global__ void PairModel_postForceKernel(Physica::PlainStruct<Core::device_obj<PairModel<Derived>>> pair) {
            pair.getDerived().postForceKernelImpl();
        }
    }

    template<class Derived>
    device_obj<PairModel<Derived>>::device_obj(size_t numParticle, ScalarType cutoff_)
            : cutoff(std::move(cutoff_))
            , cell(numParticle)
            , swapBuffer(numParticle * Dim) {
        squared_cutoff = square(cutoff_);
        if constexpr (!IsPotDependOnAtomIndex)
            pot_shift = pot_functor(0, 0, cutoff, squared_cutoff);
    }

    template<class Derived>
    __host__ __device__ inline typename device_obj<PairModel<Derived>>::ScalarType
    device_obj<PairModel<Derived>>::pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const {
        return Base::getDerived().pot_functor(i, j, r, r2);
    }

    template<class Derived>
    __host__ __device__ inline typename device_obj<PairModel<Derived>>::ScalarType
    device_obj<PairModel<Derived>>::force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const {
        return Base::getDerived().force_functor(i, j, r, r2);
    }

    template<class Derived>
    template<class Executor, bool IsSmallCell>
    Vector<typename device_obj<PairModel<Derived>>::ScalarType>
    device_obj<PairModel<Derived>>::force(const MDCellType& hostCell) {
        forceAsync<PageLockedVector, Executor, IsSmallCell>(hostCell, swapBuffer);
        StreamPool::getStream().wait();
        return swapBuffer;
    }

    template<class Derived>
    template<class VectorType, class Executor, bool IsSmallCell>
    void device_obj<PairModel<Derived>>::forceAsync(const MDCellType& hostCell, ContinuousVector<VectorType>& result) {
        static_assert(std::is_same<Executor, CudaExecutor>::value, "[Error]: Incorrect type of executor");
        StreamPool::getStream().wait(); //Ensure reentrancy
        swapBuffer = hostCell.getPos().flatten();
        auto flatten_pos = cell.getPos().flatten();
        swapBuffer.toDeviceAsync(flatten_pos);
        cell.setLattice(hostCell);

        dim3 gridDims;
        unsigned int numThread;
        if constexpr (IsSmallCell) {
            const auto& lattice = hostCell.getLattice();
            const auto range = MDCellType::estimateRange(lattice, cutoff);
            Index3D temp;
            for (int i = 0; i < Dim; ++i)
                temp[i] = static_cast<size_t>(2 * range[i] + 1);
            cellList.setCellGridDim(std::move(temp));

            constexpr size_t MaxAtomPerBlock = 128;
            const size_t numParticle = cell.getNumParticle();
            const auto numBlock = cellList.getNumCell();
            gridDims.x = (numParticle + MaxAtomPerBlock - 1) / MaxAtomPerBlock;
            gridDims.y = numBlock;
            gridDims.z = 1;
            numThread = numParticle > MaxAtomPerBlock ? MaxAtomPerBlock : numParticle;
            if (numBlock != forceBuffer.getColumn())
                forceBuffer.resize(cell.getDOF(), numBlock);
        }
        else {
            [[maybe_unused]] constexpr unsigned int UintMax = std::numeric_limits<unsigned int>::max();
            assert(cellList.getCellGridDimX() < UintMax && "[Error]: Unexpected large cell");
            assert(cellList.getCellGridDimY() < UintMax && "[Error]: Unexpected large cell");
            assert(cellList.getCellGridDimZ() < UintMax && "[Error]: Unexpected large cell");
            if (forceBuffer.getColumn() != 1)
                forceBuffer.resize(cell.getDOF(), 1);

            const CellListType hostCellList(hostCell, cutoff);
            hostCellList.toDevice(cellList);
            gridDims.x = static_cast<unsigned int>(cellList.getCellGridDimX());
            gridDims.y = static_cast<unsigned int>(cellList.getCellGridDimY());
            gridDims.z = static_cast<unsigned int>(cellList.getCellGridDimZ());
            const size_t maxNumAtomInCell = hostCellList.calcMaxNumAtomInCell();
            const size_t maxThread = Utils::DeviceProp::getInstance().getProperty(0).maxThreadsPerBlock;
            numThread = maxNumAtomInCell > maxThread ? maxThread : maxNumAtomInCell;
            assert(maxThread >= maxNumAtomInCell && "[Error]: Too many particle in the cell, performance may be pool");
        }

        Internal::PairModel_forceKernel<Derived, IsSmallCell><<<gridDims, numThread, numThread * sizeof(Vector3D), StreamPool::getStream()>>>(asStruct(*this));
        if constexpr (IsSmallCell)
            Internal::PairModel_postForceKernel<Derived><<<gridDims.x, numThread, 0, StreamPool::getStream()>>>(asStruct(*this));
        forceBuffer.col(0).toHostAsync(result);
    }

    template<class Derived>
    template<class VectorType, class Executor>
    inline void device_obj<PairModel<Derived>>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) {
        if (isSmallCell(cell))
            return forceAsync<VectorType, Executor, true>(cell, result);
        else
            return forceAsync<VectorType, Executor, false>(cell, result);
    }

    template<class Derived>
    template<class Executor>
    inline Vector<typename device_obj<PairModel<Derived>>::ScalarType>
    device_obj<PairModel<Derived>>::force_short(const MDCellType& cell) {
        if (isSmallCell(cell))
            return force<Executor, true>(cell);
        else
            return force<Executor, false>(cell);
    }

    template<class Derived>
    template<class Executor>
    inline Vector<typename device_obj<PairModel<Derived>>::ScalarType>
    device_obj<PairModel<Derived>>::force_long(const MDCellType& cell) const {
        return Vector<ScalarType>(cell.getNumParticle() * 3, 0);
    }

    template<class Derived>
    void device_obj<PairModel<Derived>>::swap(device_obj& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        cutoff.swap(obj.cutoff);
        squared_cutoff.swap(squared_cutoff);
        pot_shift.swap(pot_shift);
        cell.swap(obj.cell);
        forceBuffer.swap(obj.forceBuffer);
        swapBuffer.swap(obj.swapBuffer);
    }

    template<class Derived>
    template<bool IsSmallCell>
    __device__ void device_obj<PairModel<Derived>>::forceKernelImpl() {
        extern __shared__ Vector3D posBuffer[];
        const auto& pos = cell.getPos();
        if constexpr (IsSmallCell) {
            assert(cell.getNumParticle() < INT_MAX && "[Error]: This is not a small cell");
            const int numParticle = cell.getNumParticle();
            const int atom1 = blockIdx.x * blockDim.x + threadIdx.x;
            if (atom1 >= numParticle)
                return;

            const auto& cellGridDim = cellList.getCellGridDim();
            const Index3D cellIndex = PeriodIndex3D(blockIdx.y, cellGridDim);
            Vector3D factor{};
            for (int i = 0; i < Dim; ++i)
                factor[i] = ScalarType(int(cellIndex[i]) - int(cellGridDim[i]) / 2);
            const Vector3D delta = cell.getLattice().transpose() * factor;
            const Vector3D from = pos.row(atom1) + delta;

            const int numActiveThread = cub::min(blockDim.x, numParticle - blockIdx.x * blockDim.x);
            const int numBlock = (numParticle + numActiveThread - 1) / numActiveThread;
            Vector3D force_atom1(3, 0);
            for (int block = 0; block < numBlock; ++block) {
                const int shift = block * numActiveThread;
                __syncthreads();
                posBuffer[threadIdx.x] = pos.row(shift + threadIdx.x);
                __syncthreads();

                const int toThread = cub::min(numActiveThread, int(numParticle - block * numActiveThread));
                for (int thread = 0; thread < toThread; ++thread) {
                    const int atom2 = shift + thread;
                    const Vector3D& to = posBuffer[thread];
                    Vector3D r = to - from;
                    const ScalarType r2 = r.squaredNorm();
                    const bool isNotSelf = ScalarType(std::numeric_limits<ScalarType>::min()) < r2;
                    if (isNotSelf && r2 < squared_cutoff) {
                        const ScalarType dist = sqrt(r2);
                        const ScalarType f_norm = force_functor(atom1, atom2, dist, r2);
                        force_atom1 -= r * (f_norm / dist);
                    }
                }
            }

            for (int i = 0; i < Dim; ++i) {
                const int index = atom1 * Dim + i;
                forceBuffer(index, blockIdx.y) = force_atom1[i];
            }
        }
        else {
            const Index3D centerCell{blockIdx.x, blockIdx.y, blockIdx.z};
            const size_t cell = PeriodIndex3D(centerCell, cellList.getCellGridDim()).toIndex1D();
            const size_t numAtomInCell = cellList.getNumAtomInCell(cell);
            if (threadIdx.x >= numAtomInCell)
                return;

            const size_t atom1 = cellList.getCellAtomMap()[cellList.getCellStartOffset()[cell] + threadIdx.x];
            Vector3D force_atom1(3, 0);
            /* Current cell */ {
                const Vector3D from = pos.row(atom1);
                cellList.forAtomInCell(centerCell, [this, &pos, &force_atom1, &from, atom1](size_t atom2) {
                    auto to = pos.row(atom2);
                    Vector3D r = to.asVector() - from;
                    const ScalarType r2 = r.squaredNorm();
                    const bool isNotSelf = atom1 != atom2;
                    if (isNotSelf && r2 < squared_cutoff) {
                        const ScalarType dist = sqrt(r2);
                        const ScalarType f_norm = force_functor(atom1, atom2, dist, r2);
                        force_atom1 -= r * (f_norm / dist);
                    }
                });
            }
            cellList.forNeighInRange(centerCell, [this, &pos, &force_atom1, atom1](Vector3D translate, Index3D neigh) {
                const Vector3D from = pos.row(atom1) - translate;
                cellList.forAtomInCell(neigh, [this, &pos, &force_atom1, &from, atom1](size_t atom2) {
                    auto to = pos.row(atom2);
                    Vector3D r = to.asVector() - from;
                    const ScalarType r2 = r.squaredNorm();
                    if (r2 < squared_cutoff) {
                        const ScalarType dist = sqrt(r2);
                        const ScalarType f_norm = force_functor(atom1, atom2, dist, r2);
                        force_atom1 -= r * (f_norm / dist);
                    }
                });
            });

            for (int i = 0; i < Dim; ++i) {
                const size_t index = atom1 * Dim + i;
                forceBuffer(index, 0) = force_atom1[i];
            }
        }
    }

    template<class Derived>
    __device__ void device_obj<PairModel<Derived>>::postForceKernelImpl() {
        const unsigned int atom = blockIdx.x * blockDim.x + threadIdx.x;
        for (int i = 0; i < Dim; ++i) {
            const unsigned int index = atom * Dim + i;
            forceBuffer(index, 0) = forceBuffer.row(index).sum();
        }
    }
}
