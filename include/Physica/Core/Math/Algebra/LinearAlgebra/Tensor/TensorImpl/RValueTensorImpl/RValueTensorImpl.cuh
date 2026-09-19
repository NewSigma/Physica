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

#include "../RValueTensor.cuh"
#include "Physica/PlainStruct.h"

namespace Physica {
    template<class Derived, Scalar ScalarT>
    __host__ __device__ void device_obj<RValueTensor<Derived, ScalarT>>::assign(this const auto& self, Tensor auto&& target) {
        target.assert_assign(self);
        if (IsHost()) {
            auto func = [source_ = asStruct(self), target_ = asStruct(target)] __device__() mutable {
                const auto& source = source_.getDerived();
                auto& target = target_.getDerived();
                const size_t size = source.getSize();
                const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
                if (i < size) {
                    const auto indices = source.toIndexND(i);
                    target[indices] = source.calc(indices);
                }
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(func, self.makeKernelConfig());
        }

        if constexpr (IsDevice()) {
            const size_t size = self.getSize();
            for (size_t i = 0; i < size; ++i) {
                const auto indices = self.toIndexND(i);
                target[indices] = self.calc(indices);
            }
        }
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ void device_obj<RValueTensor<Derived, ScalarT>>::assert_assign(const Tensor auto& source) const noexcept {
        if constexpr (std::same_as<device_obj<Derived>, std::remove_cvref_t<decltype(source)>>)
            assert(this != &source && "[Error]: Self assign is likely a bug");
        assert(getShape() == source.getShape() && "[Error]: Shape mismatch between two tensors");
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ constexpr KernelConfig device_obj<RValueTensor<Derived, ScalarT>>::makeKernelConfig() const noexcept {
        const size_t size = getSize();
        const uint32_t numThread = std::min<uint32_t>(size, CUDADevAttr::DefaultThreadsPerBlock);
        const uint32_t numBlock = (size + numThread - 1) / numThread;
        return KernelConfig(numBlock, numThread);
    }

    template<class Derived, Scalar ScalarT>
    __device__ decltype(auto) device_obj<RValueTensor<Derived, ScalarT>>::calc(const IndexType& indices) const {
        return Base::getDerived().calc(indices);
    }

    template<class Derived, Scalar ScalarT>
    __device__ decltype(auto) device_obj<RValueTensor<Derived, ScalarT>>::calc(std::integral auto... dims) const {
        static_assert(sizeof...(dims) == NDim, "[Error]: NDim is not consistent");
        return calc(IndexType({static_cast<size_t>(dims)...}));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueTensor<Derived, ScalarT>>::toIndex1D(const IndexType& indices) const noexcept {
        return IndexType::toIndex1D(getShape(), indices);
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueTensor<Derived, ScalarT>>::toIndexND(size_t index) const noexcept {
        return IndexType::toIndexND(getShape(), index);
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ void device_obj<RValueTensor<Derived, ScalarT>>::resize(this auto& self, const Tensor auto& x) {
        self.resize(x.getShape());
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueTensor<Derived, ScalarT>>::resize(this auto& self, std::integral auto... dims) {
        static_assert(sizeof...(dims) == NDim, "[Error]: NDim is not consistent");
        return self.resize(IndexType({static_cast<size_t>(dims)...}));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueTensor<Derived, ScalarT>>::resize(this auto& self, IndexType shape) {
        return self.resize(std::move(shape));
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ size_t device_obj<RValueTensor<Derived, ScalarT>>::dim(int index) const noexcept {
        return Base::getDerived().dim(index);
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ auto device_obj<RValueTensor<Derived, ScalarT>>::getShape() const noexcept -> IndexType {
        return Base::getDerived().getShape();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ size_t device_obj<RValueTensor<Derived, ScalarT>>::getSize() const noexcept {
        size_t size = dim(0);
        for (int i = 1; i < NDim; ++i)
            size *= dim(i);
        return size;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueTensor<Derived, ScalarT>>::isForwardDiff() noexcept {
        return ScalarType::isForwardDiff();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueTensor<Derived, ScalarT>>::isReverseDiff() noexcept {
        return ScalarType::isReverseDiff();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueTensor<Derived, ScalarT>>::isDiffable() noexcept {
        return Diffable<T>;
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueTensor<Derived, ScalarT>>::isComplex() noexcept {
        return ScalarType::isComplex();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueTensor<Derived, ScalarT>>::isLValueTensor() noexcept {
        return Derived::isLValueTensor();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueTensor<Derived, ScalarT>>::isStrided() noexcept {
        return Derived::isStrided();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueTensor<Derived, ScalarT>>::isCompact() noexcept {
        return Derived::isCompact();
    }

    template<class Derived, Scalar ScalarT>
    __host__ __device__ consteval bool device_obj<RValueTensor<Derived, ScalarT>>::isSparse() noexcept {
        return Derived::isSparse();
    }
}
