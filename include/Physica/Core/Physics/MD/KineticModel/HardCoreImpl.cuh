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
        __global__ void binaryRun_kernel(
                ScalarType latticeSize,
                ScalarType collideFactor,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> phase_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> repMass_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> buffer_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> stepBuffer_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> velocity_,
                ScalarType deltaT) {
            constexpr unsigned int WarpSize = Physica::Utils::DeviceProp::WarpSize;
            assert(blockDim.x == WarpSize);

            const ScalarType collideStep = collideFactor * deltaT;
            auto& buffer = buffer_.getDerived();
            auto& stepBuffer = stepBuffer_.getDerived();
            auto& velocity = velocity_.getDerived();
            const size_t numParticle = buffer.getLength();

            __shared__ ScalarType trial[WarpSize + 1];
            trial[threadIdx.x + 1] = latticeSize;

            const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x - blockIdx.x;
            if (!(index < numParticle))
                return;

            auto& phase = phase_.getDerived();
            const auto& repMass = repMass_.getDerived();
            auto momentum = phase.head(numParticle);
            velocity[index] = hadamard(momentum, repMass).calc(index);

            const bool isNewRun = stepBuffer[0] == deltaT;
            const ScalarType from = isNewRun ? 0 : stepBuffer[0];
            if (isNewRun) {
                auto pos = phase.tail(numParticle);
                buffer[index] = pos[index];
            }
            ScalarType lStep = from;
            ScalarType rStep = deltaT;
            ScalarType to = deltaT;
            while (true) {
                const ScalarType step = to - from;
                trial[threadIdx.x] = (buffer + velocity * step).calc(index);
                __syncthreads();
                const bool flag = trial[threadIdx.x] > trial[threadIdx.x + 1] || (!trial[threadIdx.x].isPositive());
                const bool isCollided = __syncthreads_or(flag);
                if (isCollided)
                    rStep = to;
                else
                    lStep = to;
                to = (lStep + rStep) * 0.5;

                const bool isDone = lStep == to;
                const bool isDeltaSmallEnough = (rStep - lStep) < collideStep;
                if (isDone || isDeltaSmallEnough)
                    break;
            }
            auto pos = phase.tail(numParticle);
            pos[index] = from; // Pass from step to next kernel
            stepBuffer[blockIdx.x] = lStep;
            stepBuffer[gridDim.x + blockIdx.x] = rStep;
        }

        template<class ScalarType>
        __global__ void post_binaryRun_kernel(
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> phase_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> buffer_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> stepBuffer_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> velocity_,
                ScalarType deltaT) {
            ScalarType lStep, rStep;
            /* Find min step */ {
                auto& stepBuffer = stepBuffer_.getDerived();
                ScalarType nextFrom = stepBuffer[0];
                int index = 0;
                const int numStep = stepBuffer.getLength() / 2;
                for (int i = 1; i < numStep; ++i) {
                    if (stepBuffer[i] < nextFrom) {
                        nextFrom = stepBuffer[i];
                        index = i;
                    }
                }
                stepBuffer[0] = nextFrom;
                lStep = nextFrom;
                rStep = stepBuffer[numStep + index];
            }
            const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
            auto& buffer = buffer_.getDerived();
            const auto& velocity = velocity_.getDerived();
            const size_t numParticle = buffer.getLength();
            auto& phase = phase_.getDerived();
            auto pos = phase.tail(numParticle);
            if (index < numParticle) {
                const ScalarType from = pos[index];
                pos = buffer + velocity * (rStep - from);
                if (lStep != deltaT)
                    buffer[index] += velocity[index] * (lStep - from);
            }
        }

        template<class ScalarType>
        __global__ void handleCollision_kernel(
                ScalarType latticeSize,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> phase_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> mass_) {
            constexpr double MaxDouble = 1.797693E308;
            constexpr unsigned int WarpSize = Physica::Utils::DeviceProp::WarpSize;
            constexpr unsigned int NumThread = 2 * WarpSize;
            assert(blockDim.x == NumThread);

            auto& phase = phase_.getDerived();
            const auto& mass = mass_.getDerived();
            const size_t numParticle = phase.getLength() / 2;
            auto momentum = phase.head(numParticle);
            auto pos = phase.tail(numParticle);

            const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x - blockIdx.x;
            if (!(index < numParticle))
                return;

            __shared__ ScalarType buffer[NumThread + 1];
            buffer[threadIdx.x + 1] = MaxDouble;
            buffer[threadIdx.x] = pos[index];
            __syncthreads();
            if (buffer[threadIdx.x] > buffer[threadIdx.x + 1]) {
                const ScalarType m1 = mass[index];
                const ScalarType m2 = mass[index + 1];
                const ScalarType p1 = momentum[index];
                const ScalarType p2 = momentum[index + 1];
                const ScalarType next_p1 = ((m1 - m2) * p1 + ScalarType(2) * m1 * p2) * reciprocal(m1 + m2);
                const ScalarType next_p2 = p1 + p2 - next_p1;
                momentum[index] = next_p1;
                momentum[index + 1] = next_p2;
            }
            __syncthreads();
            if (index == 0 && !pos[0].isPositive()) {
                momentum[0] = -momentum[0];
            }
            const size_t lastParticle = numParticle - 1;
            if (index == lastParticle && pos[lastParticle] > latticeSize) {
                momentum[lastParticle] = -momentum[lastParticle];
            }
        }

        template<class ScalarType>
        __global__ void removeDrift_kernel(Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> phase_) {
            constexpr unsigned int WarpSize = Physica::Utils::DeviceProp::WarpSize;
            constexpr unsigned int NumBlock = 2 * WarpSize;
            assert(blockDim.x == NumBlock);

            auto& phase = phase_.getDerived();
            const size_t numParticle = phase.getLength() / 2;
            auto momentum = phase.head(numParticle);

            ScalarType sum = 0;
            for (int i = threadIdx.x; i < numParticle; i += NumBlock) {
                sum += momentum[i];
            }
            __shared__ ScalarType buffer[NumBlock];
            buffer[threadIdx.x] = sum;
            __syncthreads();
            for (int i = NumBlock / 2; i != 0; i /= 2) {
                if (threadIdx.x < i)
                    buffer[threadIdx.x] += buffer[threadIdx.x + i];
                __syncthreads();
            }
            const ScalarType drift = buffer[0] / ScalarType(numParticle);
            for (int i = threadIdx.x; i < numParticle; i += NumBlock) {
                momentum[i] -= drift;
            }
        }

        template<class ScalarType>
        __global__ void step_kernel(
                ScalarType latticeSize,
                ScalarType collideFactor,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> phase_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> mass_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> repMass_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> buffer_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> stepBuffer_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> velocity_,
                size_t maxHandleNum,
                ScalarType deltaT,
                size_t handleNum) {
            constexpr unsigned int WarpSize = Physica::Utils::DeviceProp::WarpSize;
            const auto& stepBuffer = stepBuffer_.getDerived();
            const size_t numParticle = mass_.getDerived().getLength();
            size_t numCollisionBlock = 0;
            for (size_t i = 0; i < stepBuffer.getLength(); ++i)
                numCollisionBlock += stepBuffer[i] != deltaT;
            const bool isDone = numCollisionBlock == 0 && handleNum != 0;
            if (isDone) {
                removeDrift_kernel<ScalarType><<<1, 64, 0, cudaStreamTailLaunch>>>(phase_);
            }
            else if (handleNum == maxHandleNum) {
                return;
            }
            else {
                {
                    const unsigned int numThread = WarpSize;
                    const unsigned int numBlock = (numParticle + WarpSize - 2) / (WarpSize - 1);
                    binaryRun_kernel<<<numBlock, numThread, 0>>>(
                            latticeSize,
                            collideFactor,
                            phase_,
                            repMass_,
                            buffer_,
                            stepBuffer_,
                            velocity_,
                            deltaT);
                }
                {
                    const unsigned int numThread = WarpSize;
                    const unsigned int numBlock = (numParticle + WarpSize - 1) / WarpSize;
                    post_binaryRun_kernel<<<numBlock, numThread, 0>>>(phase_, buffer_, stepBuffer_, velocity_, deltaT);
                }
                {
                    const unsigned int numThread = 2 * WarpSize;
                    const unsigned int numBlock = (numParticle + numThread - 2) / (numThread - 1);
                    handleCollision_kernel<<<numBlock, numThread, 0>>>(latticeSize, phase_, mass_);
                }
                step_kernel<<<1, 1, 0, cudaStreamTailLaunch>>>(
                        latticeSize,
                        collideFactor,
                        phase_,
                        mass_,
                        repMass_,
                        buffer_,
                        stepBuffer_,
                        velocity_,
                        maxHandleNum,
                        deltaT,
                        handleNum + 1);
            }
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
        constexpr unsigned int numThread = WarpSize;
        const unsigned int numBlock = (numParticle + numThread - 2) / (numThread - 1);
        stepBuffer.resize(2 * numBlock);
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
            stepBuffer = deltaT;
        }
        Internal::step_kernel<<<1, 1, 0, StreamPool::getStream()>>>(
                latticeSize,
                collideFactor,
                asStruct(d_phase),
                asStruct(mass),
                asStruct(repMass),
                asStruct(buffer),
                asStruct(stepBuffer),
                asStruct(velocity),
                maxHandleNum,
                deltaT,
                0);
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
        stepBuffer.swap(obj.stepBuffer);
        velocity.swap(obj.velocity);
        lockedBuffer.swap(obj.lockedBuffer);
        std::swap(maxHandleNum, obj.maxHandleNum);
        deltaT.swap(obj.deltaT);
    }
}
