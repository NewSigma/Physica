/*
 * Copyright 2020-2024 Weibo He.
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

#include <cassert>
#include <random>

namespace Physica::Core {
    template<class T, size_t Length, class Allocator>
    Vector<T, Length, Allocator>::Vector(size_t length_) : Storage(length_) {}

    template<class T, size_t Length, class Allocator>
    Vector<T, Length, Allocator>::Vector(size_t length_, const T& value) : Vector(length_) {
        *this = value;
    }

    template<class T, size_t Length, class Allocator>
    Vector<T, Length, Allocator>::Vector(std::initializer_list<T> list) : Storage(std::move(list)) {
        Base::makeContinuous();
    }

    template<class T, size_t Length, class Allocator>
    Vector<T, Length, Allocator>::Vector(Storage array) noexcept : Storage(std::move(array)) {
        Base::makeContinuous();
    }

    template<class T, size_t Length, class Allocator>
    template<class Derived>
    Vector<T, Length, Allocator>::Vector(const RValueVector<Derived>& v) : Storage(v.getLength()) {
        v.getDerived().assignTo(*this);
    }

    template<class T, size_t Length, class Allocator>
    Vector<T, Length, Allocator> Vector<T, Length, Allocator>::copy() const {
        if constexpr (isReverseDiff) {
            using TracerType = typename T::TracerType;
            const size_t length = Base::getLength();
            TracerType::getInstance().reserve(length);
            This result(length);
            for (size_t i = 0; i < length; ++i)
                result[i] = Base::operator[](i).copy();
            return result;
        }
        else
            return *this;
    }

    template<class T, size_t Length, class Allocator>
    Vector<T, Length, Allocator> Vector<T, Length, Allocator>::zeros(size_t len) {
        This result{};
        result.reserve(len);
        for(size_t i = 0; i < len; ++i)
            result.get_allocator().construct(result.data() + i, T(0));
        result.setLength(len);
        return result;
    }

    template<class T, size_t Length, class Allocator>
    template<class RandomGenerator>
    Vector<T, Length, Allocator> Vector<T, Length, Allocator>::random_uniform(size_t len, RandomGenerator& gen) {
        This result(len);
        result.random_uniform(gen);
        return result;
    }

    template<class T, size_t Length, class Allocator>
    template<class RandomGenerator>
    Vector<T, Length, Allocator> Vector<T, Length, Allocator>::random_uniform(
            const Vector& v1, const Vector& v2, RandomGenerator& gen) {
        assert(v1.getLength() == v2.getLength());
        This result = random_uniform(v1.getLength(), gen);
        result = v1 + hadamard((v2 - v1), result);
        return result;
    }

    template<class T, size_t Length, class Allocator>
    template<class RandomGenerator>
    Vector<T, Length, Allocator> Vector<T, Length, Allocator>::random_normal(size_t len, RandomGenerator& gen) {
        This result(len);
        result.random_normal(gen);
        return result;
    }

    template<class T, size_t Length, class Allocator>
    template<class Distribution, class RandomGenerator>
    Vector<T, Length, Allocator> Vector<T, Length, Allocator>::random_any(
            size_t len, Distribution& dist, RandomGenerator& gen) {
        This result(len);
        result.random_any(dist, gen);
        return result;
    }
    /**
     * Both \param from and \param to are included
     */
    template<class T, size_t Length, class Allocator>
    Vector<T, Length, Allocator> Vector<T, Length, Allocator>::linspace(T from, T to, size_t count) {
        assert(from < to);
        const T step = (to - from) / T(count - 1);
        Vector result = Vector(count);
        for (size_t i = 0; i < count; ++i) {
            result[i] = from;
            from += step;
        }
        return result;
    }
}
