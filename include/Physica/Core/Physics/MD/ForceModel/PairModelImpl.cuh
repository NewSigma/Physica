/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/PlainStruct.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/PeriodIndex3D.h"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "PairModel.cuh"

namespace Physica {
    namespace Internal {
        template<class T>
        __global__ void PairModel_forceKernel(Physica::PlainStruct<device_obj<T>> pair) {
            static_assert(std::is_base_of<PairModel<T>, T>::value, "[Error]: It is expected the param is a PairModel");
            pair.getDerived().forceKernelImpl();
        }

        template<class T>
        __global__ void PairModel_postForceKernel(Physica::PlainStruct<device_obj<T>> pair) {
            static_assert(std::is_base_of<PairModel<T>, T>::value, "[Error]: It is expected the param is a PairModel");
            pair.getDerived().postForceKernelImpl();
        }

        template<class T>
        __global__ void PairModel_virialKernel(Physica::PlainStruct<device_obj<T>> pair) {
            static_assert(std::is_base_of<PairModel<T>, T>::value, "[Error]: It is expected the param is a PairModel");
            pair.getDerived().virialKernelImpl();
        }

        template<class T>
        __global__ void PairModel_postVirialKernel(Physica::PlainStruct<device_obj<T>> pair) {
            static_assert(std::is_base_of<PairModel<T>, T>::value, "[Error]: It is expected the param is a PairModel");
            pair.getDerived().postVirialKernelImpl();
        }
    }

    template<class Derived>
    device_obj<PairModel<Derived>>::device_obj(size_t numParticle) : cell(numParticle), swapBuffer(numParticle * Dim) {}

    template<class Derived>
    device_obj<PairModel<Derived>>::device_obj(size_t numParticle, ScalarType cutoff_) : device_obj(numParticle) {
        setCutoff(std::move(cutoff_));
    }

    template<class Derived>
    __host__ __device__ inline auto
    device_obj<PairModel<Derived>>::pot_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const -> ScalarType {
        return Base::getDerived().pot_functor(i, j, r, r2);
    }

    template<class Derived>
    __host__ __device__ inline auto
    device_obj<PairModel<Derived>>::force_functor(size_t i, size_t j, ScalarType r, ScalarType r2) const -> ScalarType {
        return Base::getDerived().force_functor(i, j, r, r2);
    }

    template<class Derived>
    template<class Executor>
    auto device_obj<PairModel<Derived>>::force(
            const LatticeMatrix& lattice,
            const InvLatticeMatrix& invLattice,
            const PositionMatrix& cartesianPos) -> VectorND<ScalarType> {
        forceAsync<PageLockedVector, Executor>(lattice, invLattice, cartesianPos, swapBuffer);
        CUDAContext::getInstance().wait();
        return swapBuffer;
    }

    template<class Derived>
    template<class Executor>
    inline auto device_obj<PairModel<Derived>>::force(const MDCellType& hostCell) -> VectorND<ScalarType> {
        return force<Executor>(hostCell.getLattice(), hostCell.getInvLattice(), hostCell.getPos());
    }

    template<class Derived>
    template<Vector V, class Executor>
    void device_obj<PairModel<Derived>>::forceAsync(
            const LatticeMatrix& lattice,
            const InvLatticeMatrix& invLattice,
            const PositionMatrix& cartesianPos,
            ContinuousVector<V>& result) {
        static_assert(std::is_same<Executor, CUDAExecutor>::value, "[Error]: Incorrect type of executor");
        dim3 gridDims;
        size_t numThread;
        preParallel(lattice, invLattice, cartesianPos, gridDims, numThread);
        if constexpr (IsSmallCell) {
            const auto numCell = cellList.getNumCell();
            if (numCell != forceBuffer.getCol())
                forceBuffer.resize(cell.getDOF(), numCell);
        }
        else {
            if (forceBuffer.getCol() != 1)
                forceBuffer.resize(cell.getDOF(), 1);
        }

        Internal::PairModel_forceKernel<Derived>
                <<<gridDims, numThread, numThread * sizeof(DeviceVector3D), CUDAContext::getInstance()>>>(asStruct(Base::getDerived()));
        if constexpr (IsSmallCell)
            Internal::PairModel_postForceKernel<Derived><<<gridDims.x, numThread, 0, CUDAContext::getInstance()>>>(asStruct(Base::getDerived()));
        check(cudaGetLastError());
        forceBuffer.col(0).toHostAsync(result);
    }

    template<class Derived>
    template<Vector V, class Executor>
    inline void device_obj<PairModel<Derived>>::forceAsync(const MDCellType& cell, ContinuousVector<V>& result) {
        forceAsync<V, Executor>(cell.getLattice(), cell.getInvLattice(), cell.getPos(), result);
    }

    template<class Derived>
    template<class Executor>
    inline auto device_obj<PairModel<Derived>>::force_short(const MDCellType& cell) -> VectorND<ScalarType> {
        return force<Executor>(cell);
    }

    template<class Derived>
    template<class Executor>
    inline auto device_obj<PairModel<Derived>>::force_long(const MDCellType& cell) const -> VectorND<ScalarType> {
        return VectorND<ScalarType>(cell.getNumParticle() * 3, 0);
    }

    template<class Derived>
    auto device_obj<PairModel<Derived>>::virial(
            const LatticeMatrix& lattice,
            const InvLatticeMatrix& invLattice,
            const PositionMatrix& cartesianPos) -> LatticeMatrix {
        static_assert(!IsSmallCell, "[Error]: Not implemented for small cell");
        dim3 gridDims;
        size_t numThread;
        preParallel(lattice, invLattice, cartesianPos, gridDims, numThread);
        const size_t numParticle = cartesianPos.getRow();
        virialBuffer.resize(NumVirialElem, numParticle);

        Internal::PairModel_virialKernel<Derived>
                <<<gridDims, numThread, numThread * sizeof(DeviceVector3D), CUDAContext::getInstance()>>>(asStruct(Base::getDerived()));
        Internal::PairModel_postVirialKernel<Derived>
                <<<1, NumVirialElem, 0, CUDAContext::getInstance()>>>(asStruct(Base::getDerived()));
        check(cudaGetLastError());

        auto head = swapBuffer.head(NumVirialElem);
        virialBuffer.col(0).toHostAsync(head);
        CUDAContext::getInstance().wait();
        LatticeMatrix result(Dim, Dim);
        for (int r = 0; r < Dim; ++r)
            for (int c = 0; c < Dim; ++c)
                result(r, c) = swapBuffer[r * Dim + c];
        result *= reciprocal(MDCellType::getVolume(lattice)) * ScalarType(0.5);
        return result;
    }

    template<class Derived>
    inline auto device_obj<PairModel<Derived>>::virial(const MDCellType& hostCell) -> LatticeMatrix {
        return virial(hostCell.getLattice(), hostCell.getInvLattice(), hostCell.getPos());
    }

    template<class Derived>
    void device_obj<PairModel<Derived>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        cutoff.swap(obj.cutoff);
        squared_cutoff.swap(obj.squared_cutoff);
        pot_shift.swap(obj.pot_shift);
        cell.swap(obj.cell);
        forceBuffer.swap(obj.forceBuffer);
        virialBuffer.swap(obj.virialBuffer);
        swapBuffer.swap(obj.swapBuffer);
    }

    template<class Derived>
    __device__ void device_obj<PairModel<Derived>>::forceKernelImpl() {
        DeviceVector3D atomForce(3, 0);
        auto kernel = [this, &atomForce](size_t i, size_t j, DeviceVector3D r, ScalarType norm1, ScalarType norm2) {
            const ScalarType f_norm = force_functor(i, j, norm1, norm2);
            atomForce -= r * (f_norm / norm1);
        };
        const size_t atom1 = forPairInCutoff(kernel);

        for (int i = 0; i < Dim; ++i) {
            const size_t index = atom1 * Dim + i;
            if constexpr (IsSmallCell)
                forceBuffer(index, blockIdx.y) = atomForce[i];
            else
                forceBuffer(index, 0) = atomForce[i];
        }
    }

    template<class Derived>
    __device__ void device_obj<PairModel<Derived>>::postForceKernelImpl() {
        const size_t numParticle = cell.getNumParticle();
        const size_t atom = blockIdx.x * blockDim.x + threadIdx.x;
        if (atom >= numParticle)
            return;

        for (int i = 0; i < Dim; ++i) {
            const size_t index = atom * Dim + i;
            forceBuffer(index, 0) = forceBuffer.row(index).sum();
        }
    }

    template<class Derived>
    __device__ void device_obj<PairModel<Derived>>::virialKernelImpl() {
        using DeviceMatrix = DeviceMDCell::LatticeMatrix;
        DeviceMatrix atomVirial(Dim, Dim, 0);
        auto kernel = [this, &atomVirial](size_t i, size_t j, DeviceVector3D r, ScalarType norm1, ScalarType norm2) {
            const ScalarType f_norm = force_functor(i, j, norm1, norm2);
            const DeviceVector3D f = r * (f_norm / norm1);
            atomVirial += f * r.transpose();
        };
        const size_t atom1 = forPairInCutoff(kernel);
        for (int i = 0; i < Dim; ++i) {
            for (int j = 0; j < Dim; ++j) {
                const int index = i * Dim + j;
                virialBuffer(index, atom1) = atomVirial(i, j);
            }
        }
    }

    template<class Derived>
    __device__ void device_obj<PairModel<Derived>>::postVirialKernelImpl() {
        const size_t numParticle = cell.getNumParticle();
        ScalarType virial = 0;
        for (size_t i = 0; i < numParticle; ++i)
            virial += virialBuffer(threadIdx.x, i);
        virialBuffer(threadIdx.x, 0) = virial;
    }

    template<class Derived>
    void device_obj<PairModel<Derived>>::setCutoff(ScalarType cutoff_) {
        cutoff = std::move(cutoff_);
        squared_cutoff = square(cutoff);
        if constexpr (!IsPotDependOnAtomIndex) {
            constexpr int unused = 0;
            pot_shift = pot_functor(unused, unused, cutoff, squared_cutoff);
        }
    }

    template<class Derived>
    void device_obj<PairModel<Derived>>::preParallel(
            const LatticeMatrix& lattice,
            const InvLatticeMatrix& invLattice,
            const PositionMatrix& cartesianPos,
            dim3& gridDims,
            size_t& numThread) {
        CUDAContext::getInstance().wait(); //Ensure reentrancy
        swapBuffer = cartesianPos.flatten();
        auto flatten_pos = cell.getPos().flatten();
        swapBuffer.toDeviceAsync(flatten_pos);
        cell.setLattice(lattice, invLattice);
        if constexpr (IsSmallCell) {
            const auto range = MDCellType::estimateRange(lattice, cutoff);
            Index3D temp;
            for (int i = 0; i < Dim; ++i)
                temp[i] = static_cast<size_t>(2 * range[i] + 1);
            cellList.setCellGridDim(std::move(temp));

            constexpr size_t MaxAtomPerBlock = 128;
            const size_t numParticle = cell.getNumParticle();
            const size_t numBlockPerCell = (numParticle + MaxAtomPerBlock - 1) / MaxAtomPerBlock;
            const auto numCell = cellList.getNumCell();
            gridDims.x = numBlockPerCell;
            gridDims.y = numCell;
            gridDims.z = 1;
            numThread = numParticle > MaxAtomPerBlock ? MaxAtomPerBlock : numParticle;
        }
        else {
            [[maybe_unused]] constexpr unsigned int UintMax = std::numeric_limits<unsigned int>::max();
            const CellListType hostCellList(lattice, cartesianPos, cutoff);
            hostCellList.toDevice(cellList);
            assert(cellList.getCellGridDimX() < UintMax && "[Error]: Unexpected large cell");
            assert(cellList.getCellGridDimY() < UintMax && "[Error]: Unexpected large cell");
            assert(cellList.getCellGridDimZ() < UintMax && "[Error]: Unexpected large cell");
            gridDims.x = static_cast<unsigned int>(cellList.getCellGridDimX());
            gridDims.y = static_cast<unsigned int>(cellList.getCellGridDimY());
            gridDims.z = static_cast<unsigned int>(cellList.getCellGridDimZ());
            const size_t maxNumAtomInCell = hostCellList.calcMaxNumAtomInCell();
            const size_t maxThread = CUDADevAttr::MaxThreadsPerBlock;
            numThread = maxNumAtomInCell > maxThread ? maxThread : maxNumAtomInCell;
            assert(maxThread >= maxNumAtomInCell && "[Error]: Too many particle in the cell, performance may be pool");
        }
    }

    template<class Derived>
    template<class Functor>
    __device__ size_t device_obj<PairModel<Derived>>::forPairInCutoff(Functor func) const {
        extern __shared__ DeviceVector3D posBuffer[];
        const auto& pos = cell.getPos();
        if constexpr (IsSmallCell) {
            assert(cell.getNumParticle() < INT_MAX && "[Error]: This is not a small cell");
            const int numParticle = cell.getNumParticle();
            const int atom1 = blockIdx.x * blockDim.x + threadIdx.x;
            if (atom1 >= numParticle)
                asm("exit;");

            const auto& cellGridDim = cellList.getCellGridDim();
            const Index3D cellIndex = PeriodIndex3D(blockIdx.y, cellGridDim);
            DeviceVector3D factor{};
            for (int i = 0; i < Dim; ++i)
                factor[i] = ScalarType(int(cellIndex[i]) - int(cellGridDim[i]) / 2);
            const DeviceVector3D delta = cell.getLattice().transpose() * factor;
            const DeviceVector3D from = pos.row(atom1) + delta;

            const int numActiveThread = __syncthreads_count(true);
            const int numBlock = (numParticle + numActiveThread - 1) / numActiveThread;
            for (int block = 0; block < numBlock; ++block) {
                const int shift = block * numActiveThread;
                __syncthreads();
                const int index = shift + threadIdx.x;
                const bool shouldRead = index < numParticle;
                if (shouldRead)
                    posBuffer[threadIdx.x] = pos.row(index);
                const int numRead = __syncthreads_count(shouldRead);

                for (int i = 0; i < numRead; ++i) {
                    const int atom2 = shift + i;
                    const DeviceVector3D& to = posBuffer[i];
                    DeviceVector3D r = to - from;
                    const ScalarType norm2 = r.squaredNorm();
                    const bool isNotSelf = ScalarType(std::numeric_limits<ScalarType>::min()) < norm2;
                    if (isNotSelf && norm2 < squared_cutoff) {
                        const ScalarType norm1 = sqrt(norm2);
                        func(atom1, atom2, r, norm1, norm2);
                    }
                }
            }
            return atom1;
        }
        else {
            const Index3D centerCell{blockIdx.x, blockIdx.y, blockIdx.z};
            const size_t cell = PeriodIndex3D(centerCell, cellList.getCellGridDim()).toIndex1D();
            const size_t numAtomInCell = cellList.getNumAtomInCell(cell);
            if (threadIdx.x >= numAtomInCell)
                asm("exit;");

            const size_t atom1 = cellList.getCellAtomMap()[cellList.getCellStartOffset()[cell] + threadIdx.x];
            /* Current cell */ {
                const DeviceVector3D from = pos.row(atom1);
                cellList.forAtomInCell(centerCell, [this, &pos, &from, &func, atom1](size_t atom2) {
                    auto to = pos.row(atom2);
                    DeviceVector3D r = to - from;
                    const ScalarType norm2 = r.squaredNorm();
                    const bool isNotSelf = atom1 != atom2;
                    if (isNotSelf && norm2 < squared_cutoff) {
                        const ScalarType norm1 = sqrt(norm2);
                        func(atom1, atom2, r, norm1, norm2);
                    }
                });
            }
            cellList.forNeighInRange(centerCell, [this, &pos, &func, atom1](DeviceVector3D translate, Index3D neigh) {
                const DeviceVector3D from = pos.row(atom1) - translate;
                cellList.forAtomInCell(neigh, [this, &pos, &from, &func, atom1](size_t atom2) {
                    auto to = pos.row(atom2);
                    DeviceVector3D r = to - from;
                    const ScalarType norm2 = r.squaredNorm();
                    if (norm2 < squared_cutoff) {
                        const ScalarType norm1 = sqrt(norm2);
                        func(atom1, atom2, r, norm1, norm2);
                    }
                });
            });
            return atom1;
        }
    }
}
