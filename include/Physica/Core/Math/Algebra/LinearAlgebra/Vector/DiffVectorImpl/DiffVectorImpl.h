/*
 * Copyright 2024 Weibo He.
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

namespace Physica::Core {
    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(size_t length) : v(length), g(length) {}

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(size_t length, T init)
            : v(length, init), g(length, 0) {}

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(size_t length, ScalarType init) requires(isForwardDiff)
            : v(length, init.value()), g(length, init.grad()) {}

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(initializer_list list) : v(list.size()), g(list.size(), 0) {
        size_t i = 0;
        for (auto& elem : list) {
            v[i] = elem.value();
            i += 1;
        }
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(ValueVector v_, GradVector g_) noexcept
            : v(std::move(v_)), g(std::move(g_)) {}

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<Vector V>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::DenseVector(const V& v) requires(!ReverseDiff<V>) : DenseVector(v.getLength()) {
        v.assignTo(*this);
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::zero_grad() {
        g = T(0);
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    inline void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::resize(size_t size) {
        v.resize(size);
        g.resize(size);
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RandomGenerator R>
    inline void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_uniform() {
        v.template random_uniform<R>();
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RandomGenerator R>
    inline void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_normal() {
        v.template random_normal<R>();
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<class Distribution, RandomGenerator R>
    inline void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_any(Distribution& dist) {
        v.template random_any<Distribution, R>(dist);
        zero_grad();
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    void DenseVector<Diff<T, Mode, Order>, Length, Allocator>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        v.swap(obj.v);
        g.swap(obj.g);
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    __host__ __device__ inline DenseVector<Diff<T, Mode, Order>, Length, Allocator>::PtrTy
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::data_ptr(size_t index) noexcept {
        assert(index < getLength() && "[Error]: Index out of range");
        return PtrTy(v.data_ptr(index), g.data_ptr(index));
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    __host__ __device__ inline DenseVector<Diff<T, Mode, Order>, Length, Allocator>::ConstPtrTy
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::data_ptr(size_t index) const noexcept {
        return const_cast<This&>(*this).data_ptr(index);
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RandomGenerator R>
    inline DenseVector<Diff<T, Mode, Order>, Length, Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_uniform(size_t len) {
        This result(len);
        result.template random_uniform<R>();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<RandomGenerator R>
    inline DenseVector<Diff<T, Mode, Order>, Length, Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_normal(size_t len) {
        This result(len);
        result.template random_normal<R>();
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    template<class Distribution, RandomGenerator R>
    inline DenseVector<Diff<T, Mode, Order>, Length, Allocator>
    DenseVector<Diff<T, Mode, Order>, Length, Allocator>::random_any(size_t len, Distribution& dist) {
        This result(len);
        result.template random_any<Distribution, R>(dist);
        return result;
    }

    template<Scalar T, DiffMode Mode, int Order, size_t Length, class Allocator>
    auto DenseVector<Diff<T, Mode, Order>, Length, Allocator>::linspace(T from, T to, size_t count) {
        return This(ValueVector::linspace(from, to, count), GradVector(count, T(0)));
    }
}
