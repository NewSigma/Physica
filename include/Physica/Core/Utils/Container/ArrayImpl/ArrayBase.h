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

#include <cassert>
#include <cstddef>
#include <memory>
#include <type_traits>
#include "Physica/CRTPBase.h"

namespace Physica {
    /**
     * Public parts among specializations of \class Array.
     */
    template<class Derived, class Allocator>
    class ArrayBase : public CRTPBase<ArrayBase<Derived, Allocator>> {
        using This = ArrayBase<Derived, Allocator>;
        using Base = CRTPBase<This>;
        using AllocatorTraits = std::allocator_traits<Allocator>;

        template<class> class Iterator;
    public:
        using allocator_type = Allocator;
        using value_type = AllocatorTraits::value_type;
        using pointer = AllocatorTraits::pointer;
        using const_pointer = AllocatorTraits::const_pointer;
        using lvalue_reference = AllocatorTraits::lvalue_reference;
        using const_lvalue_reference = AllocatorTraits::const_lvalue_reference;
        using rvalue_reference = AllocatorTraits::rvalue_reference;

        using ElemType = Traits<Derived>::ElemType;
        static_assert(std::is_same<value_type, ElemType>::value, "[Error]: Declaration is not self consistent");
    public:
        ~ArrayBase() = default;
        /* Operators */
        [[nodiscard]] __host__ __device__ auto& operator[](this auto&&, size_t index) noexcept;
        [[nodiscard]] __host__ __device__ bool operator==(const ArrayBase& array) const;
        [[nodiscard]] __host__ __device__ bool operator!=(const ArrayBase& array) const { return !(*this == array); }
        /* Operations */
        [[nodiscard]] __host__ __device__ constexpr auto begin(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ constexpr auto end(this auto&&) noexcept;

        void send(int from, int to);
        void sendrecv(int send_to, int recv_from);
        void bcast(int root);

        void read(const auto& loc, const char* name);
        void write(auto& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return Base::getDerived().size(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return Base::getDerived().getCapacity(); }
        [[nodiscard]] __host__ __device__ bool empty() const { return getLength() == 0; }
        [[nodiscard]] __host__ __device__ bool full() const noexcept { return getLength() == getCapacity(); }
        [[nodiscard]] __host__ __device__ auto* data(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto* data_ptr(this auto&&, size_t index) noexcept;
        [[nodiscard]] __host__ __device__ auto& front(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto& back(this auto&&) noexcept;
        /* Static members */
        template<class... Args>
        [[nodiscard]] __host__ __device__ consteval static bool isTrivialDefaultConstruct() noexcept;
    protected:
        ArrayBase() = default;
        ArrayBase(const ArrayBase&) = default;
        ArrayBase(ArrayBase&&) noexcept = default;
        /* Operators */
        ArrayBase& operator=(const ArrayBase&) = default;
        ArrayBase& operator=(ArrayBase&&) noexcept = default;
    };

    template<class Derived, class Allocator>
    template<class Container>
    class ArrayBase<Derived, Allocator>::Iterator {
        using This = Iterator<Container>;
        using ElemType = Traits<std::remove_const_t<Container>>::ElemType;
        constexpr static bool isConst = std::is_const<Container>::value;
    public:
        using iterator_category = std::contiguous_iterator_tag;
        using value_type = std::conditional<isConst, const ElemType, ElemType>::type;
        using difference_type = std::ptrdiff_t;
        using pointer = std::add_pointer<value_type>::type;
        using reference = std::add_lvalue_reference<value_type>::type;
    private:
        pointer p;
    public:
        constexpr Iterator() = default;
        [[gnu::always_inline]] __host__ __device__ constexpr explicit Iterator(pointer p) noexcept;
        constexpr Iterator(const This&) noexcept = default;
        constexpr Iterator(This&&) noexcept = default;
        constexpr ~Iterator() = default;
        /* Operators */
        constexpr This& operator=(const This&) = default;
        constexpr This& operator=(This&&) noexcept = default;
        [[gnu::always_inline]] __host__ __device__ constexpr This& operator++() noexcept;
        [[gnu::always_inline]] __host__ __device__ constexpr This& operator--() noexcept;
        [[gnu::always_inline]] __host__ __device__ constexpr This& operator+=(difference_type n) noexcept;
        [[gnu::always_inline]] __host__ __device__ constexpr This& operator-=(difference_type n) noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr This operator++(int) noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr This operator--(int) noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr reference operator*() const noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr pointer operator->() const noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr reference operator[](difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr bool operator==(const This& other) const noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr auto operator<=>(const This& other) const noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr This operator+(difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr This operator-(difference_type n) const noexcept;
        [[nodiscard, gnu::always_inline]] __host__ __device__ constexpr difference_type operator-(const This& other) const noexcept;
        /* Friends */
        [[gnu::always_inline]] __host__ __device__ friend constexpr This operator+(difference_type n, const This& ite) noexcept { return ite + n; }
    };
}

namespace Physica {
    template<class T, class Allocator>
    class Traits<ArrayBase<T, Allocator>> {
    public:
        using Derived = T;
    };
}

#include "ArrayBaseImpl.h"
#ifdef PHYSICA_MPI
    #include "ArrayMPI.h"
#endif
