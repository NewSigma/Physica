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

#include "Physica/Macro.h"

namespace Physica {
    template<class Derived, class Allocator> class ArrayBase;
    /**
     * PtrIteratorF(Forward)
     */
    template<class Container>
    class PtrIteratorF : public std::contiguous_iterator_tag {
        constexpr static bool isConst = std::is_const<Container>::value;
        using This = PtrIteratorF<Container>;
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
        PtrIteratorF() = default;
        __host__ __device__ explicit PtrIteratorF(pointer p) noexcept;
        __host__ __device__ PtrIteratorF(const This& ite) noexcept;
        ~PtrIteratorF() = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
        [[nodiscard]] __host__ __device__ bool operator==(const This& other) const noexcept;
        [[nodiscard]] __host__ __device__ auto operator<=>(const This& other) const noexcept;
        [[nodiscard]] __host__ __device__ This operator+(difference_type n) const noexcept;
        __host__ __device__ This& operator+=(difference_type n) noexcept;
        [[nodiscard]] __host__ __device__ This operator-(difference_type n) const noexcept;
        __host__ __device__ This& operator-=(difference_type n) noexcept;
        [[nodiscard]] __host__ __device__ difference_type operator-(const This& ite) const noexcept;
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
     * PtrIteratorR(Reverse)
     */
    template<class Container>
    class PtrIteratorR : public std::contiguous_iterator_tag {
        constexpr static bool isConst = std::is_const<Container>::value;
        using This = PtrIteratorR<Container>;
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
        PtrIteratorR() = default;
        __host__ __device__ explicit PtrIteratorR(pointer p) noexcept;
        __host__ __device__ PtrIteratorR(const This& ite) noexcept;
        ~PtrIteratorR() = default;
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
