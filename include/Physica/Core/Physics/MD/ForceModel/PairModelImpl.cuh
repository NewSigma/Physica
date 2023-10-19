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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Grid/PeriodIndex3D.h"
#include "Physica/Core/Parallel/Executor/CudaExecutor.cuh"
#include "Physica/Utils/CUDA/DeviceProp.cuh"
#include "Physica/Utils/CUDA/PlainStruct.h"

namespace Physica::Core {
    namespace Internal {
        template<class Derived>
        __global__ void PairModel_forceKernel(Physica::PlainStruct<Core::device_obj<PairModel<Derived>>> pair) {
            pair.getDerived().forceKernelImpl();
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

        dim3 gridDims{};
        /* Make blockDim */ {
            const auto& lattice = hostCell.getLattice();
            const auto range = MDCellType::estimateRange(lattice, cutoff);
            gridDims.x = static_cast<unsigned int>(2 * range[0] + 1);
            gridDims.y = static_cast<unsigned int>(2 * range[1] + 1);
            gridDims.z = static_cast<unsigned int>(2 * range[2] + 1);
        }
        const auto numBlock = gridDims.x * gridDims.y * gridDims.z;
        if (numBlock != forceBuffer.getColumn())
            forceBuffer.resize(cell.getDOF(), numBlock);

        const size_t numParticle = cell.getNumParticle();
        const size_t maxThread = Utils::DeviceProp::getInstance().getProperty(0).maxThreadsPerBlock;
        const unsigned int numThread = numParticle > maxThread ? maxThread : numParticle;
        assert(maxThread >= numParticle && "[Error]: Too many particle in the cell, performance may be pool");
        Internal::PairModel_forceKernel<Derived><<<gridDims, numThread, 0, StreamPool::getStream()>>>(asStruct(*this));
        Internal::PairModel_postForceKernel<Derived><<<1, numThread, 0, StreamPool::getStream()>>>(asStruct(*this));
        forceBuffer.col(0).toHostAsync(result);
    }

    template<class Derived>
    void device_obj<PairModel<Derived>>::swap(device_obj& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        cutoff.swap(obj.cutoff);
        squared_cutoff.swap(squared_cutoff);
        pot_shift.swap(pot_shift);
        cell.swap(obj.cell);
        forceBuffer.swap(obj.forceBuffer);
        swapBuffer.swap(obj.swapBuffer);
    }

    template<class Derived>
    __device__ void device_obj<PairModel<Derived>>::forceKernelImpl() {
        Vector3D factor{};
        factor[0] = ScalarType(int(blockIdx.x) - int(gridDim.x) / 2);
        factor[1] = ScalarType(int(blockIdx.y) - int(gridDim.y) / 2);
        factor[2] = ScalarType(int(blockIdx.z) - int(gridDim.z) / 2);

        const Vector3D delta = cell.getLattice().transpose() * factor;
        const size_t numParticle = cell.getNumParticle();
        const auto& pos = cell.getPos();
        const unsigned int threadId = threadIdx.x;
        const Vector3D from = pos.row(threadId) + delta;

        Vector3D force_thread(3, 0);
        for (size_t j = 0; j < numParticle; ++j) {
            auto to = pos.row(j);
            Vector3D r = to.asVector() - from;
            const ScalarType r2 = r.squaredNorm();
            const bool isNotSelf = ScalarType(std::numeric_limits<ScalarType>::min()) < r2;
            if (isNotSelf && r2 < squared_cutoff) {
                const ScalarType dist = sqrt(r2);
                const ScalarType f_norm = force_functor(threadId, j, dist, r2);
                force_thread -= r * (f_norm / dist);
            }
        }

        const unsigned int blockId = (blockIdx.x * gridDim.y + blockIdx.y) * gridDim.z + blockIdx.z;
        auto bufferCol = forceBuffer.col(blockId);
        for (int i = 0; i < Dim; ++i)
            bufferCol[threadId * Dim + i] = force_thread[i];
    }

    template<class Derived>
    __device__ void device_obj<PairModel<Derived>>::postForceKernelImpl() {
        const unsigned int threadId = threadIdx.x;
        for (int i = 0; i < Dim; ++i) {
            unsigned int index = threadId * Dim + i;
            forceBuffer(index, 0) = forceBuffer.row(index).sum();
        }
    }
}
