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

#include "../LValueTensor.cuh"
#include "Physica/PlainStruct.h"

namespace Physica {
    template<class Derived>
    __host__ __device__ auto device_obj<LValueTensor<Derived>>::operator=(Scalar auto x) -> device_obj<Derived>& {
        if (IsHost()) {
            auto func = [t_ = asStruct(Base::getDerived()), x] __device__() mutable {
                auto& t = t_.getDerived();
                const size_t size = t.getSize();
                const size_t i = blockIdx.x * blockDim.x + threadIdx.x;
                if (i < size) {
                    const auto indices = t.toIndexND(i);
                    t[indices] = x;
                }
            };
            CUDAExecutor::launch<CUDADevAttr::DefaultThreadsPerBlock>(func, Base::makeKernelConfig());
        }

        if constexpr (IsDevice()) {
            auto& t = Base::getDerived();
            const size_t size = t.getSize();
            for (size_t i = 0; i < size; ++i) {
                const auto indices = t.toIndexND(i);
                t[indices] = x;
            }
        }
        return Base::getDerived();
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueTensor<Derived>>::operator+=(Scalar auto x) {
        auto& t = Base::getDerived();
        (t + x).assign(t);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueTensor<Derived>>::operator-=(Scalar auto x) {
        auto& t = Base::getDerived();
        (t - x).assign(t);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueTensor<Derived>>::operator*=(Scalar auto x) {
        auto& t = Base::getDerived();
        (t * x).assign(t);
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueTensor<Derived>>::operator/=(Scalar auto x) {
        auto& t = Base::getDerived();
        (t / x).assign(t);
    }

    template<class Derived>
    __host__ __device__ device_obj<Derived>& device_obj<LValueTensor<Derived>>::operator=(const Tensor auto& x) {
        auto& target = Base::getDerived();
        target.resize(x);
        x.assign(target);
        return target;
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueTensor<Derived>>::operator+=(const Tensor auto& x) {
        Base::getDerived() = Base::getDerived() + x;
    }

    template<class Derived>
    __host__ __device__ void device_obj<LValueTensor<Derived>>::operator-=(const Tensor auto& x) {
        Base::getDerived() += x * Trv(-1);
    }

    template<class Derived>
    __device__ decltype(auto) device_obj<LValueTensor<Derived>>::operator[](this auto&& self, const IndexType& index) {
        return *self.data_ptr(index);
    }

    template<class Derived>
    __device__ decltype(auto) device_obj<LValueTensor<Derived>>::operator[](this auto&& self, std::integral auto... dims) {
        static_assert(sizeof...(dims) == Base::NDim, "[Error]: NDim is not consistent");
        return self[IndexType({static_cast<size_t>(dims)...})];
    }

    template<class Derived>
    template<RNG R>
    void device_obj<LValueTensor<Derived>>::random_uniform() {
        Derived::template random_uniform<R>(Base::getDerived().getShape()).toDeviceAsync(Base::getDerived());
    }

    template<class Derived>
    template<RNG R>
    void device_obj<LValueTensor<Derived>>::random_normal() {
        Derived::template random_normal<R>(Base::getDerived().getShape()).toDeviceAsync(Base::getDerived());
    }

    template<class Derived>
    __device__ auto device_obj<LValueTensor<Derived>>::data_ptr(this auto&& self, const IndexType& index) noexcept {
        return self.getDerived().data_ptr(index);
    }

    template<class Derived>
    __device__ auto device_obj<LValueTensor<Derived>>::data_ptr(this auto&& self, std::integral auto... dims) noexcept {
        static_assert(sizeof...(dims) == Base::NDim, "[Error]: NDim is not consistent");
        return self.data_ptr(IndexType({static_cast<size_t>(dims)...}));
    }
}
