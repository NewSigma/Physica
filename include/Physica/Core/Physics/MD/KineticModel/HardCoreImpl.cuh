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

#include <iostream>
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

namespace Physica::Core {
    namespace Internal {
        template<class ScalarType>
        __device__ inline void handleCollision(
                typename HardCore<ScalarType, CudaExecutor>::DeviceVector& phase,
                const typename HardCore<ScalarType, CudaExecutor>::DeviceVector& mass,
                ScalarType latticeSize,
                ScalarType* sharedBuffer) {
            constexpr double MaxDouble = 1.797693E308;
            const size_t numParticle = phase.getLength() / 2;
            const unsigned int threadId = threadIdx.x;
            auto momentum = phase.head(numParticle);
            auto pos = phase.tail(numParticle);
            sharedBuffer[threadId + 1] = MaxDouble;
            sharedBuffer[threadId] = pos[threadId];
            __syncthreads();
            if (threadId < numParticle - 1 && pos[threadId] > pos[threadId + 1]) {
                const ScalarType m1 = mass[threadId];
                const ScalarType m2 = mass[threadId + 1];
                const ScalarType p1 = momentum[threadId];
                const ScalarType p2 = momentum[threadId + 1];
                const ScalarType next_p1 = ((m1 - m2) * p1 + ScalarType(2) * m1 * p2) * reciprocal(m1 + m2);
                const ScalarType next_p2 = p1 + p2 - next_p1;
                momentum[threadId] = next_p1;
                momentum[threadId + 1] = next_p2;
            }
            if (threadId == 0 && !pos[0].isPositive()) {
                momentum[0] = abs(momentum[0]);
            }
            __syncthreads();
            const size_t lastParticle = numParticle - 1;
            if (threadId && lastParticle && pos[lastParticle] > latticeSize) {
                momentum[lastParticle] = -abs(momentum[lastParticle]);
            }
            __syncthreads();
        }

        template<class ScalarType>
        __device__ inline void removeDrift(typename HardCore<ScalarType, CudaExecutor>::DeviceVector& phase, ScalarType* sharedBuffer) {
            const unsigned int numThread = blockDim.x;
            const unsigned int threadId = threadIdx.x;
            const size_t numParticle = phase.getLength() / 2;
            auto momentum = phase.head(numParticle);

            ScalarType sum = 0;
            for (int i = threadId; i < numParticle; i += numThread)
                sum += momentum[i];
            sharedBuffer[threadId] = sum;
            
            for (int shift = (numThread + 1) / 2; shift != 0; shift /= 2) {
                __syncthreads();
                if (threadId < shift && threadId + shift < numThread)
                    sharedBuffer[threadId] += sharedBuffer[threadId + shift];
            }
            __syncthreads();
            const ScalarType drift = sharedBuffer[0] / ScalarType(numParticle);
            for (int i = threadId; i < numParticle; i += numThread) {
                momentum[i] -= drift;
            }
        }

        template<class ScalarType>
        __global__ void step_kernel(
                ScalarType latticeSize,
                ScalarType collideStep,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> phase_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> mass_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> repMass_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> buffer_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> velocity_,
                ScalarType deltaT,
                size_t maxHandleNum) {
            auto& phase = phase_.getDerived();
            const auto& repMass = repMass_.getDerived();
            auto& buffer = buffer_.getDerived();
            auto& velocity = velocity_.getDerived();

            const unsigned int threadId = threadIdx.x;
            const size_t numParticle = buffer.getLength();
            if (!(threadId < numParticle))
                return;
            assert(numParticle == 512);

            auto momentum = phase.head(numParticle);
            auto pos = phase.tail(numParticle);
            __shared__ ScalarType sharedBuffer[512 + 1];
            buffer[threadId] = pos[threadId];
            velocity[threadId] = hadamard(momentum, repMass).calc(threadId);
            sharedBuffer[threadId + 1] = latticeSize;

            
            ScalarType lStep = 0;
            ScalarType rStep = deltaT;
            ScalarType from = 0;
            ScalarType to = deltaT;
            size_t handleNum = 0;
            while (true) {
                __syncthreads();
                const ScalarType step = to - from;
                sharedBuffer[threadId] = (buffer + velocity * step).calc(threadId);
                __syncthreads();
                const bool flag = sharedBuffer[threadId] > sharedBuffer[threadId + 1]
                                || (!sharedBuffer[threadId].isPositive());
                const bool isCollided = __syncthreads_or(flag);
                if (isCollided)
                    rStep = to;
                else
                    lStep = to;
                to = (lStep + rStep) * 0.5;

                const bool isDone = lStep == deltaT;
                if (isDone) {
                    pos[threadId] = sharedBuffer[threadId];
                    break;
                }
                const bool isDeltaSmallEnough = (rStep - lStep) < collideStep;
                if (isDeltaSmallEnough) {
                    if (handleNum == maxHandleNum) [[unlikely]]
                        return;
                    handleNum += 1;
                    pos[threadId] = buffer[threadId] + velocity[threadId] * (rStep - from);
                    buffer[threadId] += velocity[threadId] * (lStep - from);
                    from = lStep;
                    to = deltaT;
                    rStep = deltaT;
                    handleCollision(phase, mass_.getDerived(), latticeSize, sharedBuffer);
                    velocity[threadId] = hadamard(momentum, repMass).calc(threadId);
                    sharedBuffer[threadId + 1] = latticeSize;
                }
            }
            removeDrift(phase, sharedBuffer);
        }
    }

    template<class ScalarType>
    HardCore<ScalarType, CudaExecutor>::HardCore(ScalarType latticeSize_, ScalarType collideFactor_, size_t numParticle, size_t maxHandleNum_)
            : latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , d_phase(2 * numParticle)
            , mass(numParticle)
            , repMass(numParticle)
            , buffer(numParticle)
            , velocity(numParticle)
            , lockedBuffer(numParticle * 2)
            , maxHandleNum(maxHandleNum_)
            , deltaT(0) {
        assert(collideFactor < ScalarType(1.0) && collideFactor.isPositive());
    }

    template<class ScalarType>
    HardCore<ScalarType, CudaExecutor>& HardCore<ScalarType, CudaExecutor>::operator=(HardCore<ScalarType, CudaExecutor> obj) noexcept {
        swap(*this);
        return *this;
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::nve_step(RingPolymerType& ringPolymer, ScalarType deltaT_) {
        pre_nve_step(ringPolymer);
        do_nve_step(ringPolymer, deltaT_);
        post_nve_step(ringPolymer);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::nve_step_for(ScalarType duration, RingPolymerType& ringPolymer, ScalarType deltaT_) {
        pre_nve_step(ringPolymer);
        do_nve_step_for(duration, ringPolymer, deltaT_);
        post_nve_step(ringPolymer);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::pre_nve_step(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        lockedBuffer = phase;
        lockedBuffer.toDeviceAsync(d_phase);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::do_nve_step(RingPolymerType& ringPolymer, ScalarType deltaT_) {
        assert(getNumParticle() == ringPolymer.getNumParticle());
        assert(deltaT_.isPositive());
        if (deltaT != deltaT_) {
            deltaT = deltaT_;
        }
        const size_t numParticle = getNumParticle();
        const unsigned int numThread = numParticle > 1024 ? 1024 : numParticle;
        Internal::step_kernel<<<1, numThread, 0, StreamPool::getStream()>>>(
                latticeSize,
                collideFactor * deltaT,
                asStruct(d_phase),
                asStruct(mass),
                asStruct(repMass),
                asStruct(buffer),
                asStruct(velocity),
                deltaT,
                maxHandleNum);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::do_nve_step_for(ScalarType duration, RingPolymerType& ringPolymer, ScalarType deltaT_) {
        const uint64_t step = double(duration / deltaT_) + 0.5;
        for (uint64_t _ = 0; _ < step; ++_)
            do_nve_step(ringPolymer, deltaT_);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::post_nve_step(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        d_phase.toHostAsync(lockedBuffer);
        CudaExecutor::wait();
        phase = lockedBuffer;
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::updateMomentum(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        auto momentum = phase.head(numParticle);
        auto d_momentum = d_phase.head(numParticle);

        auto head = lockedBuffer.head(getNumParticle());
        head = momentum;
        head.toDeviceAsync(d_momentum);
        velocity = hadamard(d_momentum, repMass);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::updateMass(RingPolymerType& ringPolymer) {
        auto head = lockedBuffer.head(getNumParticle());
        head = ringPolymer.getMassVec();
        head.toDeviceAsync(mass);
        repMass = reciprocal(mass);
        CudaExecutor::wait();
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::swap(HardCore& obj) noexcept {
        latticeSize.swap(obj.latticeSize);
        collideFactor.swap(obj.collideFactor);
        d_phase.swap(obj.d_phase);
        mass.swap(obj.mass);
        repMass.swap(obj.repMass);
        buffer.swap(obj.buffer);
        velocity.swap(obj.velocity);
        lockedBuffer.swap(obj.lockedBuffer);
        std::swap(maxHandleNum, obj.maxHandleNum);
        deltaT.swap(obj.deltaT);
    }
}
