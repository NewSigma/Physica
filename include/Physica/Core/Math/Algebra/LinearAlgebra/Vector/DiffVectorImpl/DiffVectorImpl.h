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

#include "../DiffVector.h"

namespace Physica {
    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(size_t length) noexcept : v(length) {
        if constexpr (isForwardDiff())
            g.resize(length);
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(size_t length, T init) : v(length, init), g(GradVector::zeros(length)) {}

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(size_t length, ScalarType init) requires(isForwardDiff()) : v(length, init.value()), g(length, init.grad()) {}

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(std::initializer_list<T> list) : v(list) {
        if constexpr (isForwardDiff())
            g.resize(v.getLength());
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(std::initializer_list<ScalarType> list) requires(isForwardDiff()) : DenseVector(list.size()) {
        auto view = zip(v.view(), g.view(), list);
        for (auto it = view.begin(); it < view.end(); ++it) {
            auto [it_v, it_g, it_l] = it;
            *it_v = it_l->value();
            *it_g = it_l->grad();
        }
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(ValueVector v_, GradVector g_) noexcept
            : v(std::move(v_)), g(std::move(g_)) {}

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(const Vector auto& src) : DenseVector(src.getLength()) {
        using V = decltype(src);
        static_assert(!ReverseDiff<V>, "[Error]: Copying a reverse diff object breaks the compute graph");
        if constexpr (!Diffable<V>) {
            src.assign(v);
            zero_grad();
        }
        else
            src.assign(*this);
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    bool DenseVector<Diff<T, Mode, Order>, Length, Allocator>::operator==(const This& other) const {
        return v == other.v && g == other.g;
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::zero_grad() {
        g.zeros();
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::resize(size_t size) {
        v.resize(size);
        if constexpr (isForwardDiff())
            g.resize(size);
        else
            assert(g.empty() && "[Error]: Reject resize-after-reverse");
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RNG R>
    void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_uniform() {
        v.template random_uniform<R>();
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RNG R>
    void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_normal() {
        v.template random_normal<R>();
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RNG R>
    void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_any(auto& distribution) {
        v.template random_any<R>(distribution);
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        v.swap(obj.v);
        g.swap(obj.g);
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    auto DenseVector<Diff<T, Mode, Order>, Length, Allocator>::data(this auto&& self) noexcept {
        constexpr bool IsConst = std::is_const<std::remove_reference_t<decltype(self)>>::value;
        using U = Diff<T, Mode, Order>;
        using RetTy = std::conditional<IsConst, typename U::ConstPtrTy, typename U::PtrTy>::type;
        return RetTy(self.values().data(), self.grads().data());
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    auto&& DenseVector<Diff<T, Mode, Order>, Length, Allocator>::values(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.v);
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<int GradOrder>
    auto&& DenseVector<Diff<T, Mode, Order>, Length, Allocator>::grads(this auto&& self) noexcept {
        if constexpr (isReverseDiff()) {
            if (self.g.empty()) {
                auto& g1 = const_cast<GradVector&>(self.g);
                g1.resize(self.v);
                g1.zeros();
            }
        }

        if constexpr (GradOrder == 1)
            return forward_like<decltype(self)>(self.g);
        else
            return forward_like<decltype(self)>(self.g.template grads<GradOrder - 1>());
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    __host__ __device__ consteval size_t DenseVector<Diff<T, Mode, Order>, Length, Allocator>::getSizeAtCompile() noexcept {
        return Length;
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RNG R>
    auto DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_uniform(size_t len) -> This {
        This result(len);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RNG R>
    auto DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_normal(size_t len) -> This {
        This result(len);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RNG R>
    auto DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_any(size_t len, auto& distribution) -> This {
        This result(len);
        result.template random_any<R>(distribution);
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    auto DenseVector<Diff<T, Mode, Order>, Length, Allocator>::linspace(T from, T to, size_t count) {
        return This(ValueVector::linspace(from, to, count));
    }
}
