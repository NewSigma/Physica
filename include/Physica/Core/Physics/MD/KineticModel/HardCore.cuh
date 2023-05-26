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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.cuh"
#include "Physica/Core/Parallel/Executor/CudaExecutor.cuh"
#include "Physica/Core/Parallel/CudaEvent.cuh"
#include "Physica/Utils/CUDA/PlainStruct.h"
#include "Physica/Utils/CUDA/DeviceProp.cuh"
#include "Physica/Utils/Container/PageLockedAllocator.cuh"
#include "HardCore.h"

namespace Physica::Core {
    template<class ScalarType>
    class HardCore<ScalarType, CudaExecutor> {
        constexpr static unsigned int WarpSize = Physica::Utils::DeviceProp::WarpSize;
        static_assert(std::is_pointer_v<cudaGraph_t>);
    public:
        using RingPolymerType = typename HardCore<ScalarType>::RingPolymerType;
        using DeviceVector = device_obj<Vector<ScalarType>>;
        using PageLockedVector = Vector<ScalarType, Dynamic, Dynamic, Utils::PageLockedAllocator<ScalarType>>;
    private:
        ScalarType latticeSize;
        ScalarType collideFactor;
        DeviceVector d_phase;
        DeviceVector mass;
        DeviceVector repMass;
        DeviceVector buffer;
        DeviceVector stepBuffer;
        DeviceVector velocity;
        PageLockedVector lockedBuffer;
        size_t maxHandleNum;
        ScalarType deltaT;


        void* binaryRunKernelParams[8];
        void* postBinaryRunKernelParams[5];
        void* handleCollisionKernelParams[3];
        cudaKernelNodeParams binaryRunKernelNodeParams;
        cudaKernelNodeParams postBinaryRunKernelNodeParams;
        cudaKernelNodeParams handleCollisionKernelNodeParams;
        CudaEvent copyDoneEvent;

        cudaGraphExec_t binaryRunGraphExec;
        cudaGraph_t binaryRunGraph;
        cudaGraphNode_t binaryRunKernelNode;
        cudaGraphNode_t copyStepsNode;
        cudaGraphNode_t copyDoneEventNode;
        cudaGraphNode_t postBinaryRunKernelNode;
        cudaGraphNode_t handleCollisionKernelNode;

        cudaGraphExec_t binaryRunNoCopyGraphExec;
        cudaGraph_t binaryRunNoCopyGraph;
        cudaGraphNode_t binaryRunKernelNode1;
        cudaGraphNode_t postBinaryRunKernelNode1;
        cudaGraphNode_t handleCollisionKernelNode1;
    public:
        HardCore(ScalarType latticeSize_, ScalarType collideFactor_, size_t numParticle, size_t maxHandleNum_);
        HardCore(const HardCore&) = default;
        HardCore(HardCore&&) noexcept = default;
        ~HardCore();
        /* Operators */
        HardCore& operator=(HardCore obj) noexcept;
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, ScalarType deltaT_);
        void nve_step_for(ScalarType duration, RingPolymerType& ringPolymer, ScalarType deltaT_);

        void pre_nve_step(RingPolymerType& ringPolymer);
        void do_nve_step(RingPolymerType& ringPolymer, ScalarType deltaT_);
        void do_nve_step_for(ScalarType duration, RingPolymerType& ringPolymer, ScalarType deltaT_);
        void post_nve_step(RingPolymerType& ringPolymer);

        void updateMomentum(RingPolymerType& ringPolymer);
        void updateMass(RingPolymerType& ringPolymer);
        void swap(HardCore& obj) noexcept; //TODO: graph related
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getNumParticle() const noexcept { return repMass.getLength(); }
    private:
        ScalarType binaryRun();
        void handleCollision(RingPolymerType& ringPolymer);
    };
}

#include "HardCoreImpl.cuh"
