/*
 * Copyright 2020-2025 Weibo He.
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

#include "../DenseVector.h"

namespace Physica {
    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator>::DenseVector(size_t length) : Storage(length) {}

    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator>::DenseVector(size_t length, const T& init) : DenseVector(length) {
        *this = init;
    }

    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator>::DenseVector(std::initializer_list<T> list) : Storage(std::move(list)) {}

    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator>::DenseVector(Storage array) noexcept : Storage(std::move(array)) {}

    template<Scalar T, size_t Length, class Allocator>
    template<Vector V>
    DenseVector<T, Length, Allocator>::DenseVector(const V& v) : Storage(v.getLength()) {
        v.assign(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator> DenseVector<T, Length, Allocator>::zeros(size_t len) {
        This result{};
        result.reserve(len);
        for(size_t i = 0; i < len; ++i)
            result.get_allocator().construct(result.data() + i, T(0));
        result.setLength(len);
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    DenseVector<T, Length, Allocator> DenseVector<T, Length, Allocator>::random_uniform(size_t len) {
        This result(len);
        result.random_uniform<R>();
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    DenseVector<T, Length, Allocator> DenseVector<T, Length, Allocator>::random_uniform(const This& v1, const This& v2) {
        assert(v1.getLength() == v2.getLength());
        This result = random_uniform<R>(v1.getLength());
        result = v1 + hadamard((v2 - v1), result);
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    DenseVector<T, Length, Allocator> DenseVector<T, Length, Allocator>::random_normal(size_t len) {
        This result(len);
        result.random_normal<R>();
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R, class Distribution>
    DenseVector<T, Length, Allocator> DenseVector<T, Length, Allocator>::random_any(
            size_t len, Distribution& dist) {
        This result(len);
        result.random_any<R, decltype(dist)>(dist);
        return result;
    }
    /**
     * Both \param from and \param to are included
     */
    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator> DenseVector<T, Length, Allocator>::linspace(T from, T to, size_t count) {
        assert(from < to);
        const T step = (to - from) / T(count - 1);
        This result = This(count);
        for (size_t i = 0; i < count; ++i) {
            result[i] = from;
            from += step;
        }
        return result;
    }
}
