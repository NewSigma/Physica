/*
 * Copyright 2026 Weibo He.
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
