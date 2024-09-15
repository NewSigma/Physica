/*
 * Copyright 2023-2024 Weibo He.
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

#include <iostream>
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Parallel/Future/StreamFuture.cuh"
#include "Physica/Core/Physics/MD/RPMD.h"

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
        __global__ void __launch_bounds__(512, 1)
        HardCore_stepKernel(
                Physica::PlainStruct<HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>> obj,
                ScalarType deltaT,
                size_t numStep) {
            obj.getDerived().stepKernelImpl(deltaT, numStep);
        }
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::HardCore(
            ScalarType latticeSize_, ScalarType collideFactor_, ScalarType temperatureT_, size_t numParticle, size_t maxHandleNum_)
            : latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , temperatureT(temperatureT_)
            , d_phase(2 * numParticle)
            , mass(numParticle)
            , repMass(numParticle)
            , buffer(numParticle)
            , lockedBuffer(numParticle * 2)
            , maxHandleNum(maxHandleNum_) {
        HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, SequentialExecutor>::checkParam(collideFactor, 1);
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>&
    HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::operator=(HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor> obj) noexcept {
        swap(*this);
        return *this;
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::nve_step(
            RingPolymerType& ringPolymer, ScalarType deltaT) {
        pre_nve_step(ringPolymer);
        do_nve_step(deltaT, 1);
        post_nve_step(ringPolymer);
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::nve_step_for(
            ScalarType duration, RingPolymerType& ringPolymer, ScalarType deltaT) {
        pre_nve_step(ringPolymer);
        do_nve_step_for(duration, deltaT);
        post_nve_step(ringPolymer);
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::pre_nve_step(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        lockedBuffer = phase;
        lockedBuffer.toDeviceAsync(d_phase);
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::do_nve_step(ScalarType deltaT, size_t numStep) {
        assert(deltaT.isPositive());
        const size_t numParticle = getNumParticle();
        const size_t maxThread = DeviceProp::getInstance().getProperty(0).maxThreadsPerBlock;
        const unsigned int numThread = numParticle > maxThread ? maxThread : numParticle;
        Internal::HardCore_stepKernel<ScalarType, IsFixedBoundary, NumReplica>
            <<<1, numThread, numThread * sizeof(ScalarType), CUDAContext::getInstance()>>>(asStruct(*this), deltaT, numStep);
        check(cudaGetLastError());
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::do_nve_step_for(
            ScalarType duration, ScalarType deltaT) {
        const uint64_t step = RPMDBase<ScalarType>::durationToStep(duration, deltaT);
        do_nve_step(deltaT, step);
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::post_nve_step(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        d_phase.toHostAsync(lockedBuffer);
        CUDAContext::getInstance().wait();
        phase = lockedBuffer;
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::updateMomentum(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        auto momentum = phase.head(numParticle);
        auto d_momentum = d_phase.head(numParticle);

        auto head = lockedBuffer.head(getNumParticle());
        head = momentum;
        head.toDeviceAsync(d_momentum);
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::updateMass(RingPolymerType& ringPolymer) {
        auto head = lockedBuffer.head(getNumParticle());
        head = ringPolymer.getMassVec();
        head.toDeviceAsync(mass);
        Vector<ScalarType> temp = reciprocal(head);
        temp.toDeviceAsync(repMass);
        CUDAContext::getInstance().wait();
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::swap(HardCore& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        latticeSize.swap(obj.latticeSize);
        collideFactor.swap(obj.collideFactor);
        temperatureT.swap(obj.temperatureT);
        d_phase.swap(obj.d_phase);
        mass.swap(obj.mass);
        repMass.swap(obj.repMass);
        buffer.swap(obj.buffer);
        lockedBuffer.swap(obj.lockedBuffer);
        std::swap(maxHandleNum, obj.maxHandleNum);
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    __device__ inline void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::stepKernelImpl(
            ScalarType deltaT, size_t numStep) {
        const unsigned int threadId = threadIdx.x;
        const size_t numParticle = buffer.getLength();
        if (!(threadId < numParticle))
            return;
        const ScalarType collideStep = collideFactor * deltaT;
        assert(numParticle == 512);
        assert(bool(HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator>::checkStepSize(
            latticeSize, temperatureT, collideStep, mass[threadId])));

        auto momentum = d_phase.head(numParticle);
        auto pos = d_phase.tail(numParticle);
        extern __shared__ ScalarType sharedBuffer[];
        for (size_t step = 0; step < numStep; ++step) {
            buffer[threadId] = pos[threadId];
            ScalarType velocity = hadamard(momentum, repMass).calc(threadId);
            
            ScalarType lStep = 0;
            ScalarType rStep = deltaT;
            ScalarType from = 0;
            ScalarType to = deltaT;
            size_t handleNum = 0;
            while (true) {
                /* Try step */ {
                    const ScalarType stepSize = to - from;
                    const ScalarType pos1 = buffer[threadId] + velocity * stepSize;
                    sharedBuffer[threadId] = pos1;
                    __syncthreads();
                    bool flag = false;
                    if constexpr (IsFixedBoundary) {
                        if (threadId != blockDim.x - 1)
                            flag = pos1 > sharedBuffer[threadId + 1] || (!pos1.isPositive());
                        else
                            flag = pos1 > latticeSize;
                    }
                    else {
                        if (threadId != blockDim.x - 1)
                            flag = pos1 > sharedBuffer[threadId + 1];
                        else
                            flag = pos1 > sharedBuffer[0] + latticeSize;
                    }
                    const bool isCollided = __syncthreads_or(flag);
                    if (isCollided)
                        rStep = to;
                    else
                        lStep = to;
                    to = (lStep + rStep) * 0.5;
                    const bool isDone = lStep == deltaT;
                    if (isDone) {
                        pos[threadId] = pos1;
                        break;
                    }
                }
                const bool isDeltaSmallEnough = (rStep - lStep) < collideStep;
                if (isDeltaSmallEnough) {
                    if (handleNum == maxHandleNum) [[unlikely]]
                        __trap();

                    const ScalarType pos0 = buffer[threadId];
                    pos[threadId] = pos0 + velocity * (rStep - from);
                    buffer[threadId] = pos0 + velocity * (lStep - from);

                    from = lStep;
                    to = deltaT;
                    handleCollision(sharedBuffer);
                    velocity = hadamard(momentum, repMass).calc(threadId);
                    handleNum += 1;
                }
            }
        }
    }

    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    __device__ inline void HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>::handleCollision(
            __restrict__ ScalarType* sharedBuffer) {
        const unsigned int threadId = threadIdx.x;
        const size_t numParticle = mass.getLength();
        auto momentum = d_phase.head(numParticle);
        auto pos = d_phase.tail(numParticle);

        sharedBuffer[threadId] = pos[threadId];
        __syncthreads();
        const unsigned int threadId1 = (threadId + 1) % numParticle;
        bool isCollided = false;
        if (threadId != blockDim.x - 1)
            isCollided = sharedBuffer[threadId] > sharedBuffer[threadId1];
        else {
            if constexpr (!IsFixedBoundary)
                isCollided = sharedBuffer[threadId] > sharedBuffer[0] + latticeSize;
        }

        if (isCollided) {
            const ScalarType m1 = mass[threadId];
            const ScalarType m2 = mass[threadId1];
            const ScalarType p1 = momentum[threadId];
            const ScalarType p2 = momentum[threadId1];
            const ScalarType next_p1 = ((m1 - m2) * p1 + ScalarType(2) * m1 * p2) * reciprocal(m1 + m2);
            const ScalarType next_p2 = p1 + p2 - next_p1;
            momentum[threadId] = next_p1;
            momentum[threadId1] = next_p2;
        }
        if constexpr (IsFixedBoundary) {
            if (threadId == 0 && (!pos[0].isPositive()))
                momentum[0] = abs(momentum[0]);
            __syncthreads();
            const size_t lastParticle = numParticle - 1;
            if (threadId == lastParticle && (pos[lastParticle] > latticeSize)) {
                momentum[lastParticle] = -abs(momentum[lastParticle]);
            }
        }
        __syncthreads();
    }
}
