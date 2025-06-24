/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Parallel/Parallel.h"
#include "Physica/Core/Utils/CUDA/device_obj.cuh"

namespace Physica {
    /**
     * Given \tparam HostModel, \class CPUGPUModel enables collaborative computing on both CPU and GPU.
     */
    template<class HostModel, class DeviceModel>
    class CPUGPUModel {
        static_assert(!is_device_obj<HostModel>::value, "[Error]: Host model must not be device object");
        static_assert(!is_device_obj<DeviceModel>::value, "[Error]: device_obj<> is unnecessary");
        using T = HostModel::ScalarType;
        using MDCellType = HostModel::MDCellType;
        HostModel hostModel;
        Array<device_obj<DeviceModel>> deviceModels;
    public:
        template<class... Args>
        CPUGPUModel(size_t numCUDAThread, HostModel hostModel_, Args&&... deviceArgs);

        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force(const MDCellType& cell);
        template<Vector V, ExecutePolicy P>
        void forceAsync(const MDCellType& cell, ContinuousVector<V>& result);

        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_short(const MDCellType& cell) { return force<P>(cell); }
        template<ExecutePolicy P>
        [[nodiscard]] VectorND<T> force_long(const MDCellType& cell) const { return VectorND<T>(cell.getNumParticle() * 3, 0); }
        /* Getters */
        [[nodiscard]] size_t getNumCUDAThread() const noexcept { return deviceModels.getLength(); }
    };

    template<class HostModel, class DeviceModel>
    template<class... Args>
    CPUGPUModel<HostModel, DeviceModel>::CPUGPUModel(size_t numCUDAThread, HostModel hostModel_, Args&&... deviceArgs)
            : hostModel(std::move(hostModel_)) {
        assert(numCUDAThread < ThreadPool::getInstance().getNumThreads());
        deviceModels.resize(numCUDAThread, std::forward<Args>(deviceArgs)...);
    }

    template<class HostModel, class DeviceModel>
    template<ExecutePolicy P>
    auto CPUGPUModel<HostModel, DeviceModel>::force(const MDCellType& cell) -> DenseVector<T> {
        VectorND<T> result(cell.getDOF());
        forceAsync<VectorND<T>, P>(cell, result);
        CUDAContext::getInstance().wait();
        return result;
    }

    template<class HostModel, class DeviceModel>
    template<Vector V, ExecutePolicy P>
    void CPUGPUModel<HostModel, DeviceModel>::forceAsync(const MDCellType& cell, ContinuousVector<V>& result) {
        static_assert(P == GPU, "[Error]: Invalid executor");
        const auto threadId = ThreadPool::getThreadID();
        const bool useCPU = ThreadPool::isMainThread() || static_cast<size_t>(threadId) >= getNumCUDAThread();
        if (useCPU)
            result = hostModel.template force<Thread>(cell);
        else
            deviceModels[threadId].template forceAsync<V, GPU>(cell, result);
    }
}

namespace Physica {
    template<class HostModel, class DeviceModel>
    class Traits<CPUGPUModel<HostModel, DeviceModel>> : public Traits<HostModel> {};
}
