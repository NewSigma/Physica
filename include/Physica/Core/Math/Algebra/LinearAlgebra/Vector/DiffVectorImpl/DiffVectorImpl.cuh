/*
 * Copyright 2024-2026 Weibo He.
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

#include "../DiffVector.cuh"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order>
    device_obj<DenseVector<Diff<T, Mode, Order>>>::device_obj(size_t length) : v(length), g(length) {}

    template<Scalar T, DiffMode Mode, int Order>
    device_obj<DenseVector<Diff<T, Mode, Order>>>::device_obj(size_t length, T init) : device_obj(length) {
        v = init;
        g.zeros();
    }

    template<Scalar T, DiffMode Mode, int Order>
    device_obj<DenseVector<Diff<T, Mode, Order>>>::device_obj(size_t length, ScalarType init) requires(isForwardDiff)
            : v(length, init.value()), g(length, init.grad()) {}

    template<Scalar T, DiffMode Mode, int Order>
    template<Vector V>
    device_obj<DenseVector<Diff<T, Mode, Order>>>::device_obj(const V& v_) requires(!ReverseDiff<V>) : device_obj(v_.getLength()) {
        if constexpr (!Diffable<V>) {
            v_.assign(v);
            g.zeros();
        }
        else
            v_.assign(*this);
    }

    template<Scalar T, DiffMode Mode, int Order>
    void device_obj<DenseVector<Diff<T, Mode, Order>>>::zero_grad() {
        g.zeros();
    }

    template<Scalar T, DiffMode Mode, int Order>
    void device_obj<DenseVector<Diff<T, Mode, Order>>>::resize(const Vector auto& x) {
        resize(x.getLength());
    }

    template<Scalar T, DiffMode Mode, int Order>
    void device_obj<DenseVector<Diff<T, Mode, Order>>>::resize(size_t size) {
        if (size != getLength()) {
            v.resize(size);
            g.resize(size);
        }
    }

    template<Scalar T, DiffMode Mode, int Order>
    auto device_obj<DenseVector<Diff<T, Mode, Order>>>::toHost() const -> host_obj {
        host_obj result = toHostAsync();
        CUDAContext::getInstance().wait();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order>
    auto device_obj<DenseVector<Diff<T, Mode, Order>>>::toHostAsync() const -> host_obj {
        host_obj result(getLength());
        toHostAsync(result);
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order>
    void device_obj<DenseVector<Diff<T, Mode, Order>>>::toHost(host_obj& obj) const {
        toHostAsync(obj);
        CUDAContext::getInstance().wait();
    }

    template<Scalar T, DiffMode Mode, int Order>
    void device_obj<DenseVector<Diff<T, Mode, Order>>>::toHostAsync(host_obj& obj) const {
        v.toHostAsync(obj.v);
        g.toHostAsync(obj.g);
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RNG R>
    void device_obj<DenseVector<Diff<T, Mode, Order>>>::random_uniform() {
        v.template random_uniform<R>();
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RNG R>
    void device_obj<DenseVector<Diff<T, Mode, Order>>>::random_normal() {
        v.template random_normal<R>();
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RNG R>
    void device_obj<DenseVector<Diff<T, Mode, Order>>>::random_any(auto& distribution) {
        v.template random_any<R>(distribution);
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order>
    void device_obj<DenseVector<Diff<T, Mode, Order>>>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        v.swap(obj.v);
        g.swap(obj.g);
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ auto device_obj<DenseVector<Diff<T, Mode, Order>>>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.getLength() && "[Error]: Index out of range");
        constexpr bool IsConst = std::is_const<std::remove_reference_t<decltype(self)>>::value;
        using U = Diff<T, Mode, Order>;
        using RetTy = std::conditional<IsConst, typename U::ConstPtrTy, typename U::PtrTy>::type;
        return RetTy(self.v.data_ptr(index), self.g.data_ptr(index));
    }

    template<Scalar T, DiffMode Mode, int Order>
    __host__ __device__ auto&& device_obj<DenseVector<Diff<T, Mode, Order>>>::values(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.v);
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<int GradOrder>
    __host__ __device__ auto&& device_obj<DenseVector<Diff<T, Mode, Order>>>::grads(this auto&& self) noexcept {
        if constexpr (GradOrder == 1)
            return forward_like<decltype(self)>(self.g);
        else
            return forward_like<decltype(self)>(self.g.template grads<GradOrder - 1>());
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RNG R>
    auto device_obj<DenseVector<Diff<T, Mode, Order>>>::random_uniform(size_t len) -> This {
        This result(len);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RNG R>
    auto device_obj<DenseVector<Diff<T, Mode, Order>>>::random_normal(size_t len) -> This {
        This result(len);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order>
    template<RNG R>
    auto device_obj<DenseVector<Diff<T, Mode, Order>>>::random_any(size_t len, auto& distribution) -> This {
        This result(len);
        result.template random_any<R>(distribution);
        return result;
    }
}
