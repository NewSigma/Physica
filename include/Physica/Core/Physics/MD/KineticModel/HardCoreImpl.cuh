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
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> buffer_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> stepBuffer_,
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> velocity_,
                ScalarType from,
                ScalarType deltaT) {
            constexpr unsigned int WarpSize = Physica::Utils::DeviceProp::WarpSize;
            const ScalarType collideStep = collideFactor * deltaT;
            const auto& buffer = buffer_.getDerived();
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

            ScalarType lStep = from;
            ScalarType rStep = deltaT;
            ScalarType to = deltaT;
            while (true) {
                const ScalarType step = to - from;
                trial[threadIdx.x] = (buffer + velocity * step).calc(index);
                __syncwarp();
                bool isCollided;
                {
                    __shared__ bool flags[WarpSize];
                    flags[threadIdx.x] = trial[threadIdx.x] > trial[threadIdx.x + 1] || (!trial[threadIdx.x].isPositive());
                    __syncwarp();
                    for (int i = warpSize / 2; i != 0; i /= 2)
                        if (threadIdx.x < i)
                            flags[threadIdx.x] |= flags[threadIdx.x + i];
                    __syncwarp();
                    isCollided = flags[0];
                }
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
            auto& stepBuffer = stepBuffer_.getDerived();
            stepBuffer[blockIdx.x] = lStep;
            stepBuffer[gridDim.x + blockIdx.x] = rStep;
        }

        template<class ScalarType>
        __global__ void post_binaryRun_kernel(
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> phase_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> buffer_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> stepBuffer_,
                const Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> velocity_,
                ScalarType from) {
            const auto numBlock = gridDim.x;
            const auto& stepBuffer = stepBuffer_.getDerived();
            ScalarType rStep = stepBuffer[gridDim.x];
            for (int i = 1; i < numBlock; ++i)
                if (stepBuffer[gridDim.x + i] < rStep)
                    rStep = stepBuffer[gridDim.x + i];

            const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
            const auto& buffer = buffer_.getDerived();
            const auto& velocity = velocity_.getDerived();
            const size_t numParticle = buffer.getLength();
            auto& phase = phase_.getDerived();
            auto pos = phase.tail(numParticle);
            if (index < numParticle) {
                pos = buffer + velocity * (rStep - from);
            }
        }
    }

    template<class ScalarType>
    HardCore<ScalarType, CudaExecutor>::HardCore(ScalarType latticeSize_, ScalarType collideFactor_, size_t numParticle, size_t maxHandleNum_)
            : latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , d_phase(2 * numParticle)
            , repMass(numParticle)
            , buffer(numParticle)
            , stepBuffer((numParticle + WarpSize - 2) / (WarpSize - 1) * 2)
            , velocity(numParticle)
            , lockedBuffer(numParticle)
            , maxHandleNum(maxHandleNum_) {
        assert(collideFactor < ScalarType(1.0) && collideFactor.isPositive());
        cudaCheck(cudaGraphCreate(&binaryRunGraph, 0));
        const unsigned int numBlock = stepBuffer.getLength() / 2;
        {
            binaryRunKernelNodeParams.func = reinterpret_cast<void*>(Internal::binaryRun_kernel<ScalarType>);
            binaryRunKernelNodeParams.gridDim = dim3{numBlock};
            binaryRunKernelNodeParams.blockDim = dim3{WarpSize};
            binaryRunKernelNodeParams.sharedMemBytes = 0;
            binaryRunKernelParams[0] = &latticeSize;
            binaryRunKernelParams[1] = &collideFactor;
            binaryRunKernelParams[2] = &d_phase;
            binaryRunKernelParams[3] = &repMass;
            binaryRunKernelParams[4] = &buffer;
            binaryRunKernelParams[5] = &stepBuffer;
            binaryRunKernelParams[6] = &velocity;
            binaryRunKernelParams[7] = &from;
            binaryRunKernelParams[8] = &deltaT;
            binaryRunKernelNodeParams.kernelParams = &binaryRunKernelParams[0];
            binaryRunKernelNodeParams.extra = nullptr;
            cudaCheck(cudaGraphAddKernelNode(&binaryRunKernelNode, binaryRunGraph, nullptr, 0, &binaryRunKernelNodeParams));
        }
        cudaCheck(cudaGraphAddMemcpyNode1D(&copyStepsNode, binaryRunGraph, &binaryRunKernelNode, 1, lockedBuffer.data(), stepBuffer.data(), 2 * numBlock * sizeof(ScalarType), cudaMemcpyDeviceToHost));
        {
            postBinaryRunKernelNodeParams.func = reinterpret_cast<void*>(Internal::post_binaryRun_kernel<ScalarType>);
            postBinaryRunKernelNodeParams.gridDim = dim3{numBlock};
            postBinaryRunKernelNodeParams.blockDim = dim3{WarpSize};
            postBinaryRunKernelNodeParams.sharedMemBytes = 0;
            postBinaryRunKernelParams[0] = &d_phase;
            postBinaryRunKernelParams[1] = &buffer;
            postBinaryRunKernelParams[2] = &stepBuffer;
            postBinaryRunKernelParams[3] = &velocity;
            postBinaryRunKernelParams[4] = &from;
            postBinaryRunKernelNodeParams.kernelParams = &postBinaryRunKernelParams[0];
            postBinaryRunKernelNodeParams.extra = nullptr;
            cudaCheck(cudaGraphAddKernelNode(&postBinaryRunKernelNode, binaryRunGraph, &binaryRunKernelNode, 1, &postBinaryRunKernelNodeParams));
        }
        cudaCheck(cudaGraphInstantiate(&binaryRunGraphExec, binaryRunGraph, nullptr, nullptr, 0));
    }

    template<class ScalarType>
    HardCore<ScalarType, CudaExecutor>::~HardCore() {
        cudaGraphDestroy(binaryRunGraph);
        cudaGraphExecDestroy(binaryRunGraphExec);
        binaryRunGraph = nullptr;
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
        auto d_pos = d_phase.tail(numParticle);
        auto d_momentum = d_phase.head(numParticle);

        phase.toDeviceAsync(d_phase);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::do_nve_step(RingPolymerType& ringPolymer, ScalarType deltaT_) {
        const size_t numParticle = getNumParticle();
        auto d_momentum = d_phase.head(numParticle);
        auto d_pos = d_phase.tail(numParticle);
        buffer = d_pos;
        assert(numParticle == ringPolymer.getNumParticle());

        from = 0;
        deltaT = deltaT_;
        size_t handleNum = 0;
        while (true) {
            const auto lStep = binaryRun();
            const bool isDone = lStep == deltaT;
            if (isDone)
                break;
            else {
                if (handleNum == maxHandleNum) [[unlikely]]
                    throw BadConvergenceException("[Error]: Too many collision with in a step");
                handleNum += 1;
                buffer += velocity * (lStep - from);
                from = lStep;
                d_pos.toHostAsync(lockedBuffer);
                CudaExecutor::wait();
                handleCollision(ringPolymer);
                lockedBuffer.toDeviceAsync(d_momentum);
            }
        }
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
        auto pos = phase.tail(numParticle);
        auto d_pos = d_phase.tail(numParticle);
        d_pos.toHostAsync(pos);
        CudaExecutor::wait();
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::updateMomentum(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        auto momentum = phase.head(numParticle);
        auto d_momentum = d_phase.head(numParticle);

        momentum.toDeviceAsync(d_momentum);
        velocity = hadamard(d_momentum, repMass);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::updateMass(RingPolymerType& ringPolymer) {
        ringPolymer.getMassVec().toDeviceAsync(repMass);
        repMass = reciprocal(repMass);
        CudaExecutor::wait();
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::swap(HardCore& obj) noexcept {
        latticeSize.swap(obj.latticeSize);
        collideFactor.swap(obj.collideFactor);
        d_phase.swap(obj.d_phase);
        repMass.swap(obj.repMass);
        buffer.swap(obj.buffer);
        stepBuffer.swap(obj.stepBuffer);
        velocity.swap(obj.velocity);
        lockedBuffer.swap(obj.lockedBuffer);
        std::swap(maxHandleNum, obj.maxHandleNum);
        std::swap(binaryRunGraphExec, obj.binaryRunGraphExec);
        std::swap(binaryRunGraph, obj.binaryRunGraph);
        std::swap(binaryRunKernelParams, obj.binaryRunKernelParams);
        std::swap(binaryRunKernelNodeParams, obj.binaryRunKernelNodeParams);
        std::swap(binaryRunKernelNode, obj.binaryRunKernelNode);
        std::swap(copyStepsNode, obj.copyStepsNode);
        std::swap(postBinaryRunKernelParams, obj.postBinaryRunKernelParams);
        std::swap(postBinaryRunKernelNodeParams, obj.postBinaryRunKernelNodeParams);
        std::swap(postBinaryRunKernelNode, obj.postBinaryRunKernelNode);
    }

    template<class ScalarType>
    ScalarType HardCore<ScalarType, CudaExecutor>::binaryRun() {
        const size_t numParticle = getNumParticle();
        int device;
        cudaGetDevice(&device);
        const int numThread = WarpSize;
        const int numBlock = stepBuffer.getLength() / 2;
        cudaCheck(cudaGraphExecKernelNodeSetParams(binaryRunGraphExec, binaryRunKernelNode, &binaryRunKernelNodeParams));
        cudaCheck(cudaGraphExecKernelNodeSetParams(binaryRunGraphExec, postBinaryRunKernelNode, &postBinaryRunKernelNodeParams));
        cudaCheck(cudaGraphLaunch(binaryRunGraphExec, StreamPool::getStream()));

        CudaExecutor::wait();

        auto temp1 = lockedBuffer.head(2 * numBlock);
        const auto lStep = temp1.head(numBlock).min();
        return lStep;
    }
    /**
     * Special implementation to make full use of limited pinned memory.
     */
    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::handleCollision(RingPolymerType& ringPolymer) {
        const size_t numParticle = ringPolymer.getNumParticle();
        auto phase = ringPolymer.asMatrix().col(0);
        auto pos = phase.tail(numParticle);
        bool isDrifted = false;
        if (!lockedBuffer[0].isPositive()) {
            phase[0].toOpposite();
            isDrifted = true;
        }

        const auto& mass = ringPolymer.getMassVec();
        size_t i = 0;
        for (; i < numParticle - 1; ++i) {
            if (lockedBuffer[i] > lockedBuffer[i + 1]) {
                const ScalarType m1 = mass[i];
                const ScalarType m2 = mass[i + 1];
                const ScalarType p1 = phase[i];
                const ScalarType p2 = phase[i + 1];
                const ScalarType next_p1 = ((m1 - m2) * p1 + ScalarType(2) * m1 * p2) * reciprocal(m1 + m2);
                const ScalarType next_p2 = p1 + p2 - next_p1;
                phase[i] = next_p1;
                phase[i + 1] = next_p2;
            }
        }

        if (lockedBuffer[i] > latticeSize) {
            phase[i].toOpposite();
            isDrifted = true;
        }
        if (isDrifted)
            ringPolymer.removeDrift();
        lockedBuffer = phase.head(numParticle);
    }
}
