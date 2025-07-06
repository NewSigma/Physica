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

#include <cstddef>
#include <type_traits>
#include <memory>
#include "Physica/CRTPBase.h"
#include "Iterator.h"

namespace Physica {
    /**
     * Public parts among specializations of \class Array.
     */
    template<class Derived, class Allocator>
    class ArrayBase : public CRTPBase<ArrayBase<Derived, Allocator>> {
        using This = ArrayBase<Derived, Allocator>;
        using Base = CRTPBase<This>;
        using AllocatorTraits = std::allocator_traits<Allocator>;
    public:
        using allocator_type = Allocator;
        using value_type = AllocatorTraits::value_type;
        using pointer = AllocatorTraits::pointer;
        using const_pointer = AllocatorTraits::const_pointer;
        using lvalue_reference = AllocatorTraits::lvalue_reference;
        using const_lvalue_reference = AllocatorTraits::const_lvalue_reference;
        using rvalue_reference = AllocatorTraits::rvalue_reference;

        using ElemType = Traits<Derived>::ElemType;
    protected:
        using FIteType = FIterator<Derived>;
        using RIteType = RIterator<Derived>;
        using CFIteType = FIterator<const Derived>;
        using CRIteType = RIterator<const Derived>;

        static_assert(std::is_same<value_type, ElemType>::value, "[Error]: Declaration is not self consistent");
    public:
        ~ArrayBase() = default;
        /* Operators */
        [[nodiscard]] __host__ __device__ lvalue_reference operator[](size_t index);
        [[nodiscard]] __host__ __device__ const_lvalue_reference operator[](size_t index) const;
        __host__ __device__ bool operator==(const ArrayBase& array) const;
        __host__ __device__ bool operator!=(const ArrayBase& array) const { return !(*this == array); }
        /* Iterators */
        __host__ __device__ FIteType begin() noexcept { return FIteType(data()); }
        __host__ __device__ CFIteType begin() const noexcept { return cbegin(); }
        __host__ __device__ CFIteType cbegin() const noexcept { return CFIteType(data()); }
        __host__ __device__ FIteType end() noexcept { return FIteType(data() + getLength()); }
        __host__ __device__ CFIteType end() const noexcept { return cend(); }
        __host__ __device__ CFIteType cend() const noexcept { return CFIteType(data() + getLength()); }
        __host__ __device__ RIteType rbegin() noexcept { return RIteType(data() + getLength() - 1); }
        __host__ __device__ CRIteType rbegin() const noexcept { return crbegin(); }
        __host__ __device__ CRIteType crbegin() const noexcept { return CRIteType(data() + getLength() - 1); }
        __host__ __device__ RIteType rend() noexcept { return RIteType(data() - 1); }
        __host__ __device__ CRIteType rend() const noexcept { return crend(); }
        __host__ __device__ CRIteType crend() const noexcept { return CRIteType(data() - 1); }
        /* Operations */
        void read(const auto& loc, const char* name);
        void write(auto& loc, const char* name) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t size() const noexcept { return Base::getDerived().size(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Base::getDerived().getLength(); }
        [[nodiscard]] __host__ __device__ size_t getCapacity() const noexcept { return Base::getDerived().getCapacity(); }
        [[nodiscard]] __host__ __device__ bool empty() const { return getLength() == 0; }
        [[nodiscard]] __host__ __device__ bool full() const noexcept { return getLength() == getCapacity(); }
        [[nodiscard]] __host__ __device__ pointer data() noexcept { return Base::getDerived().data(); }
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept { return Base::getDerived().data(); }
        [[nodiscard]] __host__ __device__ pointer data_ptr(size_t index) noexcept;
        [[nodiscard]] __host__ __device__ const_pointer data_ptr(size_t index) const noexcept;
        /* Static members */
        template<class... Args>
        __host__ __device__ consteval static bool isTrivialDefaultConstruct();
    protected:
        ArrayBase() = default;
        ArrayBase(const ArrayBase&) = default;
        ArrayBase(ArrayBase&&) noexcept = default;
        /* Operators */
        ArrayBase& operator=(const ArrayBase&) = default;
        ArrayBase& operator=(ArrayBase&&) noexcept = default;
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
