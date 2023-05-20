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
#include "Physica/Utils/CUDA/PlainStruct.h"
#include "Physica/Utils/CUDA/DeviceProp.cuh"
#include "HardCore.h"

namespace Physica::Core {
    template<class ScalarType>
    class HardCore<ScalarType, CudaExecutor> {
    public:
        using RingPolymerType = typename HardCore<ScalarType>::RingPolymerType;
        using DeviceVector = device_obj<Vector<ScalarType>>;
    private:
        ScalarType latticeSize;
        ScalarType collideFactor;
        DeviceVector d_phase;
        DeviceVector repMass;
        DeviceVector buffer;
        DeviceVector velocity;
    public:
        HardCore(ScalarType latticeSize_, ScalarType collideFactor_, size_t numParticle);
        HardCore(const HardCore&) = default;
        HardCore(HardCore&&) noexcept = default;
        ~HardCore() = default;
        /* Operators */
        HardCore& operator=(HardCore obj) noexcept;
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, ScalarType deltaT);
        void updateMass(RingPolymerType& ringPolymer);
        void swap(HardCore& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumParticle() const noexcept { return repMass.getLength(); }
    private:
        bool checkCollision() const;
    };

    template<class ScalarType>
    HardCore<ScalarType, CudaExecutor>::HardCore(ScalarType latticeSize_, ScalarType collideFactor_, size_t numParticle)
            : latticeSize(latticeSize_)
            , collideFactor(collideFactor_)
            , d_phase(2 * numParticle)
            , repMass(numParticle)
            , buffer(numParticle)
            , velocity(numParticle) {
        assert(collideFactor < ScalarType(1.0) && collideFactor.isPositive());
    }

    template<class ScalarType>
    HardCore<ScalarType, CudaExecutor>& HardCore<ScalarType, CudaExecutor>::operator=(HardCore<ScalarType, CudaExecutor> obj) noexcept {
        swap(*this);
        return *this;
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::nve_step(RingPolymerType& ringPolymer, ScalarType deltaT) {
        const size_t numParticle = ringPolymer.getNumParticle();
        const ScalarType collideStep = collideFactor * deltaT;
        auto phase = ringPolymer.asMatrix().col(0);
        auto momentum = phase.head(numParticle);
        auto pos = phase.tail(numParticle);
        auto d_momentum = d_phase.head(numParticle);
        auto d_pos = d_phase.tail(numParticle);
        assert(numParticle == repMass.getLength());

        ScalarType lStep = 0;
        ScalarType rStep = deltaT;
        ScalarType from = 0;
        ScalarType to = deltaT;
        phase.toDevice(d_phase);
        buffer = d_pos;
        velocity = hadamard(d_momentum, repMass);
        while (lStep != deltaT) {
            const bool isDeltaSmallEnough = (rStep - lStep) < collideStep;
            if (isDeltaSmallEnough) {
                d_pos = buffer + velocity * (rStep - from);
                buffer += velocity * (lStep - from);
                from = lStep;
                to = deltaT;
                rStep = deltaT;
                d_pos.toHost(pos);
                HardCore<ScalarType>::handleCollision(latticeSize, ringPolymer);
                momentum.toDevice(d_momentum);
                velocity = hadamard(d_momentum, repMass);
            }

            const ScalarType step = to - from;
            d_pos = buffer + velocity * step;
            if (checkCollision())
                rStep = to;
            else
                lStep = to;
            to = (lStep + rStep) * 0.5;
        }
        d_pos.toHost(pos);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::updateMass(RingPolymerType& ringPolymer) {
        ringPolymer.getMassVec().toDevice(repMass);
        repMass = reciprocal(repMass);
    }

    template<class ScalarType>
    void HardCore<ScalarType, CudaExecutor>::swap(HardCore& obj) noexcept {
        latticeSize.swap(obj.latticeSize);
        collideFactor.swap(obj.collideFactor);
        d_phase.swap(obj.d_phase);
        repMass.swap(obj.repMass);
        buffer.swap(obj.buffer);
        velocity.swap(obj.velocity);
    }

    namespace Internal {
        template<class ScalarType>
        __global__ void checkCollision_kernel(
                Physica::PlainStruct<typename HardCore<ScalarType, CudaExecutor>::DeviceVector> phase,
                ScalarType latticeSize,
                size_t numParticle) {
            constexpr unsigned int WarpSize = Physica::Utils::DeviceProp::WarpSize;
            const unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
            auto pos = phase.getDerived().tail(numParticle);

            __shared__ ScalarType buffer[WarpSize + 1];
            buffer[index + 1] = latticeSize;
            bool flag = false;
            {
                size_t i = 0;
                const size_t delta = warpSize - 1;
                const size_t to = numParticle - delta;
                for (; i < to; i += delta) {
                    buffer[index] = pos[i + index];
                    __syncwarp();
                    flag |= buffer[index] > buffer[index + 1];
                }

                buffer[index] = latticeSize;
                if (i + index < numParticle) {
                    buffer[index] = pos[i + index];
                    __syncwarp();
                    flag |= buffer[index] > buffer[index + 1];
                }
            }
            __shared__ bool flags[WarpSize];
            flags[index] = flag;
            __syncwarp();
            for (int i = warpSize / 2; i != 0; i /= 2)
                if (index < i)
                    flags[index] |= flags[index + i];
            if (index == 0 && flags[0])
                pos[0] = -latticeSize;
        }
    }

    template<class ScalarType>
    bool HardCore<ScalarType, CudaExecutor>::checkCollision() const {
        using namespace Physica;
        int device;
        cudaGetDevice(&device);
        const int numBlock = 1;
        const int numThread = Utils::DeviceProp::getInstance().getProperty(device).warpSize;
        Internal::checkCollision_kernel<<<numBlock, numThread, 0, StreamPool::getStream()>>>(asStruct(d_phase), latticeSize, getNumParticle());
        Vector<ScalarType, 1> temp{};
        d_phase.segment(getNumParticle(), getNumParticle() + 1).toHost(temp);
        return temp[0].isNegative();
    }
}
