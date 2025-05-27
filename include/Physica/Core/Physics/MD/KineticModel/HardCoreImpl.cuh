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
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.cuh"
#include "Physica/Core/Physics/MD/RPMD.h"

namespace Physica {
    namespace Internal {
        template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
        __global__ void __launch_bounds__(512, 1)
        HardCore_stepKernel(
                Physica::PlainStruct<HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>> obj,
                T deltaT,
                size_t numStep) {
            obj.getDerived().stepKernelImpl(deltaT, numStep);
        }
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::HardCore(
            T latticeSize_, T collideFactor_, T temperatureT_, size_t numParticle, size_t maxHandleNum_)
            : latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , temperatureT(temperatureT_)
            , d_phase(2 * numParticle)
            , mass(numParticle)
            , repMass(numParticle)
            , buffer(numParticle)
            , lockedBuffer(numParticle * 2)
            , maxHandleNum(maxHandleNum_) {
        HardCore<T, IsFixedBoundary, NumReplica, Integrator, Sequential>::checkParam(collideFactor, 1);
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::nve_step(
            RingPolymerType& ringPolymer, T deltaT) {
        pre_nve_step(ringPolymer);
        do_nve_step(deltaT, 1);
        post_nve_step(ringPolymer);
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::nve_step_for(
            T duration, RingPolymerType& ringPolymer, T deltaT) {
        pre_nve_step(ringPolymer);
        do_nve_step_for(duration, deltaT);
        post_nve_step(ringPolymer);
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::pre_nve_step(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        lockedBuffer = phase;
        lockedBuffer.toDeviceAsync(d_phase);
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::do_nve_step(T deltaT, size_t numStep) {
        assert(deltaT.isPositive());
        const size_t numParticle = getNumParticle();
        const size_t maxThread = CUDADevAttr::MaxThreadsPerBlock;
        const unsigned int numThread = numParticle > maxThread ? maxThread : numParticle;
        Internal::HardCore_stepKernel<T, IsFixedBoundary, NumReplica>
            <<<1, numThread, numThread * sizeof(T), CUDAContext::getInstance()>>>(asStruct(*this), deltaT, numStep);
        check(cudaGetLastError());
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::do_nve_step_for(
            T duration, T deltaT) {
        const uint64_t step = RPMD<T>::durationToStep(duration, deltaT);
        do_nve_step(deltaT, step);
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::post_nve_step(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        d_phase.toHostAsync(lockedBuffer);
        CUDAContext::getInstance().wait();
        phase = lockedBuffer;
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::updateMomentum(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        auto momentum = phase.head(numParticle);
        auto d_momentum = d_phase.head(numParticle);

        auto head = lockedBuffer.head(getNumParticle());
        head = momentum;
        head.toDeviceAsync(d_momentum);
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::updateMass(RingPolymerType& ringPolymer) {
        auto head = lockedBuffer.head(getNumParticle());
        head = ringPolymer.getMassVec();
        head.toDeviceAsync(mass);
        VectorND<T> temp = reciprocal(head);
        temp.toDeviceAsync(repMass);
        CUDAContext::getInstance().wait();
    }

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::swap(This& __restrict obj) noexcept {
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

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    __device__ inline void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::stepKernelImpl(
            T deltaT, size_t numStep) {
        const unsigned int threadId = threadIdx.x;
        const size_t numParticle = buffer.getLength();
        if (!(threadId < numParticle))
            return;
        const T collideStep = collideFactor * deltaT;
        assert(numParticle == 512);
        assert(bool(HardCore<T, IsFixedBoundary, NumReplica, Integrator>::checkStepSize(
            latticeSize, temperatureT, collideStep, mass[threadId])));

        auto momentum = d_phase.head(numParticle);
        auto pos = d_phase.tail(numParticle);
        extern __shared__ T sharedBuffer[];
        for (size_t step = 0; step < numStep; ++step) {
            buffer[threadId] = pos[threadId];
            T velocity = hadamard(momentum, repMass).calc(threadId);
            
            T lStep = 0;
            T rStep = deltaT;
            T from = 0;
            T to = deltaT;
            size_t handleNum = 0;
            while (true) {
                /* Try step */ {
                    const T stepSize = to - from;
                    const T pos1 = buffer[threadId] + velocity * stepSize;
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

                    const T pos0 = buffer[threadId];
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

    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    __device__ inline void HardCore<T, IsFixedBoundary, NumReplica, Integrator, GPU>::handleCollision(T* __restrict sharedBuffer) {
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
            const T m1 = mass[threadId];
            const T m2 = mass[threadId1];
            const T p1 = momentum[threadId];
            const T p2 = momentum[threadId1];
            const T next_p1 = ((m1 - m2) * p1 + T(2) * m1 * p2) * reciprocal(m1 + m2);
            const T next_p2 = p1 + p2 - next_p1;
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
