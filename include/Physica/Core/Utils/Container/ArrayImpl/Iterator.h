/*
 * Copyright 2024-2025 Weibo He.
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

#include "Physica/Macro.h"

namespace Physica {
    template<class Derived, class Allocator> class ArrayBase;
    /**
     * FIterator(Forward)
     */
    template<class Container>
    class FIterator : public std::contiguous_iterator_tag {
        constexpr static bool isConst = std::is_const<Container>::value;
        using This = FIterator<Container>;
        using ElemType = Traits<std::remove_const_t<Container>>::ElemType;
    public:
        using difference_type = std::ptrdiff_t;
        using value_type = std::conditional<isConst, const ElemType, ElemType>::type;
        using pointer = std::add_pointer<value_type>::type;
        using reference = std::add_lvalue_reference<value_type>::type;
        using iterator_category = std::contiguous_iterator_tag;
    private:
        pointer p;
    public:
        FIterator() = default;
        __host__ __device__ explicit FIterator(pointer p) noexcept;
        __host__ __device__ FIterator(const This& ite) noexcept;
        ~FIterator() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        [[nodiscard]] __host__ __device__ This operator+(difference_type n) const noexcept;
        [[nodiscard]] __host__ __device__ This& operator+=(difference_type n) const noexcept;
        [[nodiscard]] __host__ __device__ This operator-(difference_type n) const noexcept;
        [[nodiscard]] __host__ __device__ This& operator-=(difference_type n) const noexcept;
        [[nodiscard]] __host__ __device__ difference_type operator-(const This& ite) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator<(const This& ite) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator>(const This& ite) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator<=(const This& ite) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator>=(const This& ite) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator==(const This& ite) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator!=(const This& ite) const noexcept;
        __host__ __device__ This& operator++() noexcept;
        __host__ __device__ This operator++(int) noexcept;
        __host__ __device__ This& operator--() noexcept;
        __host__ __device__ This operator--(int) noexcept;
        [[nodiscard]] __host__ __device__ reference operator*() const noexcept;
        [[nodiscard]] __host__ __device__ pointer operator->() const noexcept;
        [[nodiscard]] __host__ __device__ reference operator[](difference_type n) const noexcept;
        /* Friends */
        friend __host__ __device__ This operator+(difference_type n, const This& it) noexcept { return it + n; }
    };
    /**
     * RIterator(Reverse)
     */
    template<class Container>
    class RIterator : public std::contiguous_iterator_tag {
        constexpr static bool isConst = std::is_const<Container>::value;
        using This = RIterator<Container>;
        using ElemType = Traits<std::remove_const_t<Container>>::ElemType;
    public:
        using difference_type = std::ptrdiff_t;
        using value_type = std::conditional<isConst, const ElemType, ElemType>::type;
        using pointer = std::add_pointer<value_type>::type;
        using reference = std::add_lvalue_reference<value_type>::type;
        using iterator_category = std::contiguous_iterator_tag;
    private:
        pointer p;
    public:
        RIterator() = default;
        __host__ __device__ explicit RIterator(pointer p) noexcept;
        __host__ __device__ RIterator(const This& ite) noexcept;
        ~RIterator() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        [[nodiscard]] __host__ __device__ bool operator==(const This& ite) const noexcept;
        [[nodiscard]] __host__ __device__ bool operator!=(const This& ite) const noexcept;
        __host__ __device__ This& operator++() noexcept;
        __host__ __device__ This operator++(int) noexcept;
        [[nodiscard]] __host__ __device__ reference operator*() const noexcept;
    };
}

#include "IteratorImpl.h"
