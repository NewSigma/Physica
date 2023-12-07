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
#include "Physica/Utils/Allocator/PageLockedAllocator.cuh"
#include "HardCore.h"

namespace Physica::Core {
    template<class ScalarType, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    class HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator, CudaExecutor> {
        constexpr static unsigned int WarpSize = Physica::Utils::DeviceProp::WarpSize;
        static_assert(NumReplica == 1, "[Error]: PIMD is not implemented");
        static_assert(Integrator == RPMDIntegrator::Exact, "[Error]: Cayley integrator not implemented");
    public:
        using RingPolymerType = typename HardCore<ScalarType, IsFixedBoundary, NumReplica, Integrator>::RingPolymerType;
        using DeviceVector = device_obj<Vector<ScalarType>>;
        using PageLockedVector = Vector<ScalarType, Dynamic, Dynamic, Utils::PageLockedAllocator<ScalarType>>;
    private:
        ScalarType latticeSize;
        ScalarType collideFactor;
        ScalarType temperatureT;
        DeviceVector d_phase;
        DeviceVector mass;
        DeviceVector repMass;
        DeviceVector buffer;
        DeviceVector stepBuffer;
        PageLockedVector lockedBuffer;
        size_t maxHandleNum;
    public:
        HardCore(ScalarType latticeSize_, ScalarType collideFactor_, ScalarType temperatureT_, size_t numParticle, size_t maxHandleNum_);
        HardCore(const HardCore&) = default;
        HardCore(HardCore&&) noexcept = default;
        ~HardCore() = default;
        /* Operators */
        HardCore& operator=(HardCore obj) noexcept;
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, ScalarType deltaT);
        void nve_step_for(ScalarType duration, RingPolymerType& ringPolymer, ScalarType deltaT);

        void pre_nve_step(RingPolymerType& ringPolymer);
        void do_nve_step(ScalarType deltaT, size_t numStep);
        void do_nve_step_for(ScalarType duration, ScalarType deltaT);
        void post_nve_step(RingPolymerType& ringPolymer);

        void updateMomentum(RingPolymerType& ringPolymer);
        void updateMass(RingPolymerType& ringPolymer);
        void swap(HardCore& __restrict obj) noexcept;

        __device__ inline void stepKernelImpl(ScalarType deltaT, size_t numStep);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getNumParticle() const noexcept { return repMass.getLength(); }
    private:
        __device__ inline void handleCollision(__restrict__ ScalarType* sharedBuffer);
    };
}

#include "HardCoreImpl.cuh"
