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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.cuh"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Physica/Core/Utils/Allocator/PageLockedAllocator.cuh"
#include "HardCore.h"

namespace Physica {
    template<Scalar T, bool IsFixedBoundary, size_t NumReplica, RPMDIntegrator Integrator>
    class HardCore<T, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor> {
        using This = HardCore<T, IsFixedBoundary, NumReplica, Integrator, CUDAExecutor>;
        static_assert(NumReplica == 1, "[Error]: PIMD is not implemented");
        static_assert(Integrator == RPMDIntegrator::Exact, "[Error]: Cayley integrator not implemented");
        constexpr static unsigned int WarpSize = Physica::CUDADevAttr::WarpSize;
    public:
        using RingPolymerType = HardCore<T, IsFixedBoundary, NumReplica, Integrator>::RingPolymerType;
        using DeviceVector = device_obj<VectorND<T>>;
        using PageLockedVector = DenseVector<T, Dynamic, PageLockedAllocator<T>>;
    private:
        T latticeSize;
        T collideFactor;
        T temperatureT;
        DeviceVector d_phase;
        DeviceVector mass;
        DeviceVector repMass;
        DeviceVector buffer;
        DeviceVector stepBuffer;
        PageLockedVector lockedBuffer;
        size_t maxHandleNum;
    public:
        HardCore(T latticeSize_, T collideFactor_, T temperatureT_, size_t numParticle, size_t maxHandleNum_);
        HardCore(const This&) = default;
        HardCore(This&&) noexcept = default;
        ~HardCore() = default;
        /* Operators */
        This& operator=(HardCore obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void nve_step(RingPolymerType& ringPolymer, T deltaT);
        void nve_step_for(T duration, RingPolymerType& ringPolymer, T deltaT);

        void pre_nve_step(RingPolymerType& ringPolymer);
        void do_nve_step(T deltaT, size_t numStep);
        void do_nve_step_for(T duration, T deltaT);
        void post_nve_step(RingPolymerType& ringPolymer);

        void updateMomentum(RingPolymerType& ringPolymer);
        void updateMass(RingPolymerType& ringPolymer);
        void swap(This& __restrict obj) noexcept;

        __device__ inline void stepKernelImpl(T deltaT, size_t numStep);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getNumParticle() const noexcept { return repMass.getLength(); }
    private:
        __device__ inline void handleCollision(T* __restrict sharedBuffer);
    };
}

#include "HardCoreImpl.cuh"
