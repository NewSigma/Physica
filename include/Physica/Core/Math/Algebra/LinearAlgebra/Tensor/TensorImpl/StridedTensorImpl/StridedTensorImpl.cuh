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

#include "../StridedTensor.cuh"

namespace Physica {
    template<class Derived>
    __host__ __device__ void device_obj<StridedTensor<Derived>>::assign(this const auto& self, Tensor auto&& target) {
        using Self = std::remove_cvref_t<decltype(self)>;
        using Target = std::remove_cvref_t<decltype(target)>;
        target.assert_assign(self);
        // TODO: Lower this and zeros to a for loop of strided-matrix assign/zeros once a reshape/flatten-like tensor-to-matrix view (torch.reshape style) is available
        if constexpr (host_obj::template isMinorCompact<Self, Target>()) {
            const int minorDim = std::max(host_obj::getCompactMinorDim(self), host_obj::getCompactMinorDim(target));
            const size_t size = self.getSize();
            const size_t stride = std::min(host_obj::getCompactMinorSize(self), host_obj::getCompactMinorSize(target));
            for (size_t i = 0; i < size; i += stride) {
                const auto indices = self.toIndexND(i);
                auto src = host_obj::makeCompactBlock(self, indices, minorDim);
                auto tgt = host_obj::makeCompactBlock(target, indices, minorDim);
                src.flatten().assign(tgt.flatten());
            }
        }
        else
            self.device_obj<LValueTensor<Derived>>::assign(target);
    }

    template<class Derived>
    __host__ __device__ void device_obj<StridedTensor<Derived>>::zeros() noexcept {
        if constexpr (host_obj::template isMinorCompact<Derived, Derived>()) {
            auto& self = Base::getDerived();
            const int minorDim = host_obj::getCompactMinorDim(self);
            const size_t size = self.getSize();
            const size_t stride = host_obj::getCompactMinorSize(self);
            for (size_t i = 0; i < size; i += stride) {
                const auto indices = self.toIndexND(i);
                host_obj::makeCompactBlock(self, indices, minorDim).flatten().zeros();
            }
        }
        else
            Base::zeros();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedTensor<Derived>>::getStrides() const noexcept {
        return Base::getDerived().getStrides();
    }

    template<class Derived>
    __host__ __device__ size_t device_obj<StridedTensor<Derived>>::getStride(size_t dim) const noexcept {
        return getStrides()[dim];
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedTensor<Derived>>::data_handle() noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedTensor<Derived>>::data_handle() const noexcept {
        return Base::getDerived().data_handle();
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedTensor<Derived>>::data_ptr(this auto&& self, const IndexType& index) noexcept {
        const IndexType strides = self.getStrides();
        size_t offset = 0;
        for (int i = 0; i < Base::NDim; ++i) {
            assert(index[i] < self.dim(i));
            offset += index[i] * strides[i];
        }
        return self.data_handle() + offset;
    }

    template<class Derived>
    __host__ __device__ auto device_obj<StridedTensor<Derived>>::data_ptr(this auto&& self, std::same_as<size_t> auto... dims) noexcept {
        static_assert(sizeof...(dims) == Base::NDim, "[Error]: NDim is not consistent");
        return self.data_ptr(IndexType({dims...}));
    }

    template<class Derived>
    __host__ __device__ consteval auto device_obj<StridedTensor<Derived>>::getStrideAtCompile() noexcept -> IndexType {
        return Derived::getStrideAtCompile();
    }

    template<class Derived>
    __host__ __device__ consteval auto device_obj<StridedTensor<Derived>>::getStrideAtCompile(size_t dim) noexcept -> size_t {
        return getStrideAtCompile()[dim];
    }
}
