/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../CompactTensor.cuh"
#include "Physica/Core/Exception/CUDA/CUDA.cuh"
#include "Physica/Core/Parallel/CUDAContext.cuh"

namespace Physica {
    template<class Derived>
    __host__ __device__ void device_obj<CompactTensor<Derived>>::assign(this const auto& self, Tensor auto&& target) {
        using Self = std::remove_cvref_t<decltype(self)>;
        using Target = std::remove_cvref_t<decltype(target)>;
        constexpr bool FastCopy = Self::isCompact() && Target::isCompact()
                               && !Self::isDiffable() && !Target::isDiffable()
                               && std::same_as<typename Self::ScalarType, typename Target::ScalarType>
                               && std::is_trivially_copyable_v<typename Self::ScalarType>;
        target.assert_assign(self);
        if constexpr (FastCopy && IsHost()) {
            const size_t size = self.getSize() * sizeof(typename Self::ScalarType);
            check(cudaMemcpyAsync(target.data(), self.data(), size, cudaMemcpyKind::cudaMemcpyDeviceToDevice, CUDAContext::getInstance()));
        }
        else
            self.device_obj<StridedTensor<Derived>>::assign(target);
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactTensor<Derived>>::data(this auto&& self) noexcept {
        return self.data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactTensor<Derived>>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<CompactTensor<Derived>>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ constexpr auto device_obj<CompactTensor<Derived>>::getStrides() const noexcept {
        return host_obj::getStridesImpl(Base::getDerived().getShape());
    }
}
