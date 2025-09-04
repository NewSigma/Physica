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
    DenseVector<T, Length, Allocator>::DenseVector(const Vector auto& v) : Storage(v.getLength()) {
        v.assign(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    void DenseVector<T, Length, Allocator>::resize(const Vector auto& x) {
        resize(x.getLength());
    }

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::zeros(size_t len) -> This {
        This result(len);
        result.zeros();
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    auto DenseVector<T, Length, Allocator>::random_uniform(size_t len) -> This {
        This result(len);
        result.random_uniform<R>();
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    auto DenseVector<T, Length, Allocator>::random_uniform(const This& v1, const This& v2) -> This {
        assert(v1.getLength() == v2.getLength());
        This result = random_uniform<R>(v1.getLength());
        result = v1 + hadamard((v2 - v1), result);
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    auto DenseVector<T, Length, Allocator>::random_normal(size_t len) -> This {
        This result(len);
        result.random_normal<R>();
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    template<RNG R>
    auto DenseVector<T, Length, Allocator>::random_any(size_t len, auto& distribution) -> This {
        This result(len);
        result.random_any<R>(distribution);
        return result;
    }
    /**
     * Both \param from and \param to are included
     */
    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::linspace(T from, T to, size_t count) -> This {
        assert(from < to);
        assert(count > 1);
        Trv n = Trv(count - 1);
        This result = This(count);
        for (size_t i = 0; i < count; ++i) {
            Trv factor = Trv(i) / n;
            result[i] = from * (Trv(1) - factor) + to * factor;
        }
        return result;
    }
    /**
     * Helper function that communicates with C libraries.
     */
    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::read(size_t length, const T* __restrict p) -> This {
        return This(Storage::read(length, p));
    }
}
