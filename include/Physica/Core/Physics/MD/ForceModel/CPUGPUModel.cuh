/*
 * Copyright 2024 WeiBo He.
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

#include "Physica/Core/Parallel/ThreadPool.h"
#include "Physica/Core/Parallel/Executor/AutoExecutor.h"
#include "Physica/Utils/CUDA/device_obj.cuh"

namespace Physica::Core {
    /**
     * Given \tparam HostModel, \class CPUGPUModel enables collaborative computing on both CPU and GPU.
     */
    template<class HostModel, class DeviceModel>
    class CPUGPUModel {
        static_assert(!Utils::is_device_obj<HostModel>::value, "[Error]: Host model must not be device object");
        static_assert(!Utils::is_device_obj<DeviceModel>::value, "[Error]: device_obj<> is unnecessary");
        using ScalarType = typename HostModel::ScalarType;
        using MDCellType = typename HostModel::MDCellType;
        HostModel hostModel;
        Physica::Utils::Array<device_obj<DeviceModel>> deviceModels;
    public:
        template<class... Args>
        CPUGPUModel(size_t numCudaThread, HostModel hostModel_, Args&&... deviceArgs);

        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell);
        template<class VectorType, class Executor>
        void forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result);

        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) { return force<Executor>(cell); }
        template<class Executor>
        [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getNumParticle() * 3, 0); }
        /* Getters */
        [[nodiscard]] size_t getNumCudaThread() const noexcept { return deviceModels.getLength(); }
    };

    template<class HostModel, class DeviceModel>
    template<class... Args>
    CPUGPUModel<HostModel, DeviceModel>::CPUGPUModel(size_t numCudaThread, HostModel hostModel_, Args&&... deviceArgs)
            : hostModel(std::move(hostModel_)) {
        assert(numCudaThread < ThreadPool::getInstance().getNumThreads());
        deviceModels.resize(numCudaThread, std::forward<Args>(deviceArgs)...);
    }

    template<class HostModel, class DeviceModel>
    template<class Executor>
    Vector<typename CPUGPUModel<HostModel, DeviceModel>::ScalarType> CPUGPUModel<HostModel, DeviceModel>::force(const MDCellType& cell) {
        Vector<ScalarType> result(cell.getDOF());
        forceAsync<Vector<ScalarType>, Executor>(cell, result);
        StreamPool::getStream().wait();
        return result;
    }

    template<class HostModel, class DeviceModel>
    template<class VectorType, class Executor>
    void CPUGPUModel<HostModel, DeviceModel>::forceAsync(const MDCellType& cell, ContinuousVector<VectorType>& result) {
        static_assert(Traits<Executor>::isCudaEnabled, "[Error]: Invalid executor");
        const auto threadId = ThreadPool::getThreadInfo().id;
        const bool useCPU = ThreadPool::isMainThread() || threadId >= getNumCudaThread();
        if (useCPU)
            result = hostModel.template force<ThreadExecutor>(cell);
        else
            deviceModels[threadId].template forceAsync<VectorType, CudaExecutor>(cell, result);
    }
}

namespace Physica {
    template<class HostModel, class DeviceModel>
    class Traits<Core::CPUGPUModel<HostModel, DeviceModel>> : public Traits<HostModel> {};
}
