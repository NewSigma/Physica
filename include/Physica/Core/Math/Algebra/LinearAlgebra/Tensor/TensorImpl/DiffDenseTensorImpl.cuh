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

#include "../DiffDenseTensor.cuh"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    __host__ __device__ device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::device_obj(IndexType shape) : v(shape), g(std::move(shape)) {}

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    __host__ __device__ device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::device_obj(IndexType shape, T init)
            : v(std::move(shape), init), g(v.getShape(), T(0)) {}

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    __host__ __device__ device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::device_obj(ValueTensor v_, GradTensorType g_)
            : v(std::move(v_)), g(std::move(g_)) {}

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<Tensor X>
    device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::device_obj(const X& x) requires(!ReverseDiff<X>)
            : device_obj(x.getShape()) {
        if constexpr (!Diffable<X>) {
            x.assign(v);
            g.zeros();
        }
        else
            x.assign(*this);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    __host__ __device__ void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::resize(this auto& self, IndexType shape) {
        self.v.resize(shape);
        self.g.resize(std::move(shape));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::reserve(this auto& self, size_t size) noexcept {
        self.v.reserve(size);
        self.g.reserve(size);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::toHost() const -> host_obj {
        host_obj result = toHostAsync();
        CUDAContext::getInstance().wait();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::toHostAsync() const -> host_obj {
        host_obj result(getShape());
        toHostAsync(result);
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::toHost(host_obj& obj) const {
        toHostAsync(obj);
        CUDAContext::getInstance().wait();
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::toHostAsync(host_obj& obj) const {
        v.toHostAsync(obj.values());
        g.toHostAsync(obj.grads());
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::random_uniform(this auto& self) {
        self = random_uniform<R>(self.getShape());
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::random_normal(this auto& self) {
        self = random_normal<R>(self.getShape());
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::random_any(this auto& self, auto& distribution) {
        self = random_any<R>(self.getShape(), distribution);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::zero_grad(this auto& self) {
        self.g.zeros();
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::junk(this auto& self) {
        self.v.junk();
        self.g.junk();
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        v.swap(obj.v);
        g.swap(obj.g);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    __host__ __device__ auto device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::data_handle(this auto&& self) noexcept {
        constexpr bool IsConst = std::is_const<std::remove_reference_t<decltype(self)>>::value;
        using U = Diff<T, Mode, Order>;
        using RetTy = std::conditional<IsConst, typename U::ConstPtrTy, typename U::PtrTy>::type;
        return RetTy(self.v.data(), self.g.data());
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    __host__ __device__ auto&& device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::values(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.v);
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<int GradOrder>
    __host__ __device__ auto&& device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::grads(this auto&& self) noexcept {
        if constexpr (GradOrder == 1)
            return forward_like<decltype(self)>(self.g);
        else
            return forward_like<decltype(self)>(self.g.template grads<GradOrder - 1>());
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    auto device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::random_uniform(IndexType shape) -> This {
        return This(ValueTensor::template random_uniform<R>(std::move(shape)));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    auto device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::random_normal(IndexType shape) -> This {
        return This(ValueTensor::template random_normal<R>(std::move(shape)));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    template<RNG R>
    auto device_obj<DenseTensor<Diff<T, Mode, Order>, Dims...>>::random_any(IndexType shape, auto& distribution) -> This {
        return This(ValueTensor::template random_any<R>(std::move(shape), distribution));
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::toDevice() const {
        auto result = toDeviceAsync();
        CUDAExecutor::wait();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    auto DenseTensor<Diff<T, Mode, Order>, Dims...>::toDeviceAsync() const {
        device_obj<This> result{};
        toDeviceAsync(result);
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::toDevice(device_obj<This>& obj) const {
        toDeviceAsync(obj);
        CUDAExecutor::wait();
    }

    template<Scalar T, DiffMode Mode, int Order, int... Dims>
    void DenseTensor<Diff<T, Mode, Order>, Dims...>::toDeviceAsync(device_obj<This>& obj) const {
        v.toDeviceAsync(obj.values());
        g.toDeviceAsync(obj.grads());
    }
}
