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

namespace Physica::Core {
    template<class Derived, class Allocator> class ArrayBase;
    /**
     * FIterator(Forward)
     */
    template<class Container>
    class FIterator {
        constexpr static bool isConst = std::is_const<Container>::value;
        using This = FIterator<Container>;
        using ElemType = typename Traits<std::remove_const_t<Container>>::ElemType;
    public:
        using difference_type = std::ptrdiff_t;
        using value_type = typename std::conditional<isConst, const ElemType, ElemType>::type;
        using pointer = typename std::add_pointer<value_type>::type;;
        using reference = typename std::add_lvalue_reference<value_type>::type;
        using iterator_category = std::random_access_iterator_tag;
    private:
        pointer p;
    public:
        __host__ __device__ explicit FIterator(pointer p) : p(p) {}
        __host__ __device__ FIterator(const This& ite) : p(ite.p) {}
        ~FIterator() = default;
        /* Operators */
        This& operator=(const This&) = default;
        [[nodiscard]] __host__ __device__ This operator+(difference_type n) const { return This(p + n); }
        [[nodiscard]] __host__ __device__ This operator-(difference_type n) const { return This(p - n); }
        [[nodiscard]] __host__ __device__ bool operator<(const This& ite) const { return p < ite.p; }
        [[nodiscard]] __host__ __device__ bool operator>(const This& ite) const { return p > ite.p; }
        [[nodiscard]] __host__ __device__ bool operator<=(const This& ite) const { return !(p > ite.p); }
        [[nodiscard]] __host__ __device__ bool operator>=(const This& ite) const { return !(p < ite.p); }
        [[nodiscard]] __host__ __device__ difference_type operator-(const This& ite) const { return p - ite.p; }
        [[nodiscard]] __host__ __device__ bool operator==(const This& ite) const noexcept { return p == ite.p; }
        [[nodiscard]] __host__ __device__ bool operator!=(const This& ite) const noexcept { return p != ite.p; }
        __host__ __device__ This& operator++();
        __host__ __device__ const This operator++(int);
        __host__ __device__ This& operator--();
        [[nodiscard]] __host__ __device__ reference operator*() const { return *p; }
    };
    /**
     * RIterator(Reverse)
     */
    template<class Container>
    class RIterator {
        constexpr static bool isConst = std::is_const<Container>::value;
        using This = RIterator<Container>;
        using ElemType = typename Traits<std::remove_const_t<Container>>::ElemType;
    public:
        using difference_type = std::ptrdiff_t;
        using value_type = typename std::conditional<isConst, const ElemType, ElemType>::type;
        using pointer = typename std::add_pointer<value_type>::type;;
        using reference = typename std::add_lvalue_reference<value_type>::type;
        using iterator_category = std::random_access_iterator_tag;
    private:
        pointer p;
    public:
        __host__ __device__ explicit RIterator(pointer p) : p(p) {}
        __host__ __device__ RIterator(const This& ite) : p(ite.p) {}
        ~RIterator() = default;
        /* Operators */
        This& operator=(const This&) = default;
        [[nodiscard]] __host__ __device__ bool operator==(const This& ite) const noexcept { return p == ite.p; }
        [[nodiscard]] __host__ __device__ bool operator!=(const This& ite) const noexcept { return p != ite.p; }
        __host__ __device__ This& operator++();
        __host__ __device__ const This operator++(int);
        [[nodiscard]] __host__ __device__ reference operator*() const { return *p; }
    };
}

#include "IteratorImpl.h"
