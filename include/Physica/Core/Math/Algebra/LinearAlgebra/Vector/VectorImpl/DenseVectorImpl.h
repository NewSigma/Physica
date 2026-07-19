/*
 * Copyright 2020-2026 Weibo He.
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
    DenseVector<T, Length, Allocator>::DenseVector(Array<T, Length, Allocator> storage) noexcept : storage(std::move(storage)) {}

    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator>::DenseVector(size_t length) noexcept : storage(length) {}

    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator>::DenseVector(size_t length, const T& init) : DenseVector(length) {
        *this = init;
    }

    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator>::DenseVector(std::initializer_list<T> list) : storage(list) {}

    template<Scalar T, size_t Length, class Allocator>
    DenseVector<T, Length, Allocator>::DenseVector(const Vector auto& v) : storage(v.getLength()) {
        v.assign(*this);
    }

    template<Scalar T, size_t Length, class Allocator>
    bool DenseVector<T, Length, Allocator>::operator==(const This& other) const noexcept {
        return storage == other.storage;
    }

    template<Scalar T, size_t Length, class Allocator>
    template<size_t I>
    constexpr auto&& DenseVector<T, Length, Allocator>::get(this auto&& self) noexcept {
        return forward_like<decltype(self)>(self.storage).template get<I>();
    }

    template<Scalar T, size_t Length, class Allocator>
    constexpr auto DenseVector<T, Length, Allocator>::begin(this auto&& self) noexcept {
        return self.storage.begin();
    }

    template<Scalar T, size_t Length, class Allocator>
    constexpr auto DenseVector<T, Length, Allocator>::end(this auto&& self) noexcept {
        return self.storage.end();
    }

    template<Scalar T, size_t Length, class Allocator>
    void DenseVector<T, Length, Allocator>::append(auto&&... args) noexcept {
        storage.append(std::forward<decltype(args)>(args)...);
    }

    template<Scalar T, size_t Length, class Allocator>
    void DenseVector<T, Length, Allocator>::resize(const Vector auto& x) {
        resize(x.getLength());
    }

    template<Scalar T, size_t Length, class Allocator>
    void DenseVector<T, Length, Allocator>::resize(size_t size, auto&&... args) noexcept {
        storage.resize(size, std::forward<decltype(args)>(args)...);
    }

    template<Scalar T, size_t Length, class Allocator>
    void DenseVector<T, Length, Allocator>::reserve(size_t size) noexcept {
        storage.reserve(size);
    }

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::skew() const -> T {
        This temp = *this;
        temp.standardize();
        temp = hadamard(square(temp), temp);
        const size_t length = getLength();
        const T factor = T(length * length) / T((length - 1) * (length - 2));
        return temp.mean() * factor;
    }

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::excess_kurt() const -> T {
        This temp = *this;
        temp.standardize();
        temp = square(temp);
        const T mean2 = temp.mean();
        temp = square(temp);
        const T mean1 = temp.mean();

        const size_t length = getLength();
        const T factor2 = T(length * length * 3) / T((length - 2) * (length - 3));
        const T factor1 = T(length * length * (length + 1)) / T((length - 1) * (length - 2) * (length - 3));
        return factor1 * mean1 - factor2 * mean2;
    }

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::kurt() const -> T {
        return excess_kurt() + Trv(3);
    }
    /**
     * Both \param from and \param to are included
     */
    template<Scalar T, size_t Length, class Allocator>
    void DenseVector<T, Length, Allocator>::linspace(T from, T to) {
        const size_t length = getLength();
        assert(length > 1);
        Trv n = Trv(length - 1);
        for (size_t i = 0; i < length; ++i) {
            Trv factor = Trv(i) / n;
            (*this)[i] = from * (Trv(1) - factor) + to * factor;
        }
    }

    template<Scalar T, size_t Length, class Allocator>
    void DenseVector<T, Length, Allocator>::zeros() noexcept {
        storage.zeros();
    }

    template<Scalar T, size_t Length, class Allocator>
    void DenseVector<T, Length, Allocator>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        storage.swap(obj.storage);
    }

    template<Scalar T, size_t Length, class Allocator>
    __host__ __device__ consteval size_t DenseVector<T, Length, Allocator>::getSizeAtCompile() noexcept {
        return Length;
    }

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::zeros(size_t len) -> This {
        This result(len);
        result.zeros();
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    auto* DenseVector<T, Length, Allocator>::data(this auto&& self) noexcept {
        return self.storage.data();
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

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::linspace(T from, T to, size_t count) -> This {
        This result(count);
        result.linspace(from, to);
        return result;
    }

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::read_hdf5(const H5Loc& loc, const char* name) -> This {
        This result{};
        result.read(loc, name);
        return result;
    }
    /**
     * Helper function that communicates with C libraries.
     */
    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::read(size_t length, const T* __restrict p) noexcept -> This {
        return This(Storage::read(length, p));
    }

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::generate(std::invocable<size_t> auto fn) -> This {
        static_assert(Length != Dynamic, "[Error]: Cannot infer length");
        return This(Storage::generate(std::move(fn)));
    }

    template<Scalar T, size_t Length, class Allocator>
    auto DenseVector<T, Length, Allocator>::generate(std::invocable<size_t> auto fn, size_t length) -> This {
        return This(Storage::generate(std::move(fn), length));
    }
}
