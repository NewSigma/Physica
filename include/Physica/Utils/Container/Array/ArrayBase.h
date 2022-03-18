/*
 * Copyright 2020-2021 WeiBo He.
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
#include <utility>
#include <memory>
#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Utils::Internal {
    //Forward declaration
    template<class T> class Traits;
    template<class Derived, class Allocator> class ArrayBase;

    template<class T, class Container>
    class ContainerIterator;
    
    template<class ValueType, class Derived, class Allocator>
    class ContainerIterator<ValueType, ArrayBase<Derived, Allocator>> {
    public:
        using difference_type = std::ptrdiff_t;
        using value_type = ValueType;
        using pointer = typename std::add_pointer<ValueType>::type;;
        using reference = typename std::add_lvalue_reference<ValueType>::type;
        using iterator_category = std::random_access_iterator_tag;
    private:
        pointer p;
    public:
        __host__ __device__ ContainerIterator(const ContainerIterator& ite) : p(ite.p) {}
        ~ContainerIterator() = default;
        /* Operators */
        __host__ __device__ ContainerIterator& operator=(const ContainerIterator& ite);
        __host__ __device__ difference_type operator-(const ContainerIterator& ite) { return p - ite.p; }
        __host__ __device__ [[nodiscard]] bool operator==(const ContainerIterator& ite) const noexcept { return p == ite.p; }
        __host__ __device__ [[nodiscard]] bool operator!=(const ContainerIterator& ite) const noexcept { return p != ite.p; }
        __host__ __device__ ContainerIterator& operator++();
        __host__ __device__ const ContainerIterator operator++(int);
        __host__ __device__ [[nodiscard]] reference operator*() const { return *p; }
    private:
        __host__ __device__ explicit ContainerIterator(pointer p) : p(p) {}

        friend class ArrayBase<Derived, Allocator>;
    };

    template<class Pointer, class Container>
    class ReverseContainerIterator;
    
    template<class ValueType, class Derived, class Allocator>
    class ReverseContainerIterator<ValueType, ArrayBase<Derived, Allocator>> {
    public:
        using difference_type = std::ptrdiff_t;
        using value_type = ValueType;
        using pointer = typename std::add_pointer<ValueType>::type;;
        using reference = typename std::add_lvalue_reference<ValueType>::type;
        using iterator_category = std::random_access_iterator_tag;
    private:
        pointer p;
    public:
        __host__ __device__ ReverseContainerIterator(const ReverseContainerIterator& ite) : p(ite.p) {}
        ~ReverseContainerIterator() = default;
        /* Operators */
        __host__ __device__ ReverseContainerIterator& operator=(const ReverseContainerIterator& ite);
        __host__ __device__ [[nodiscard]] bool operator==(const ReverseContainerIterator& ite) const noexcept { return p == ite.p; }
        __host__ __device__ [[nodiscard]] bool operator!=(const ReverseContainerIterator& ite) const noexcept { return p != ite.p; }
        __host__ __device__ ReverseContainerIterator& operator++();
        __host__ __device__ const ReverseContainerIterator operator++(int);
        __host__ __device__ [[nodiscard]] reference operator*() const { return *p; }
    private:
        __host__ __device__ explicit ReverseContainerIterator(pointer p) : p(p) {}

        friend class ArrayBase<Derived, Allocator>;
    };
    /**
     * Public parts among specializations of \class Array.
     */
    template<class Derived, class Allocator>
    class ArrayBase : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using allocator_type = Allocator;
        using AllocatorTraits = std::allocator_traits<allocator_type>;
        using ValueType = typename AllocatorTraits::value_type;
        using PointerType = typename AllocatorTraits::pointer;
        using const_pointer = typename AllocatorTraits::const_pointer;
        using LValueReferenceType = typename AllocatorTraits::lvalue_reference;
        using ConstLValueReferenceType = typename AllocatorTraits::const_lvalue_reference;
        using RValueReferenceType = typename AllocatorTraits::rvalue_reference;
        using Iterator = ContainerIterator<ValueType, ArrayBase<Derived, Allocator>>;
        using ConstIterator = ContainerIterator<const ValueType, ArrayBase<Derived, Allocator>>;
        using ReverseIterator = ReverseContainerIterator<ValueType, ArrayBase<Derived, Allocator>>;
        using ConstReverseIterator = ReverseContainerIterator<const ValueType, ArrayBase<Derived, Allocator>>;
    public:
        /* Operators */
        [[nodiscard]] __host__ __device__ LValueReferenceType operator[](size_t index);
        [[nodiscard]] __host__ __device__ ConstLValueReferenceType operator[](size_t index) const;
        bool operator==(const ArrayBase& array) const;
        __host__ __device__ bool operator!=(const ArrayBase& array) const { return !(*this == array); }
        /* Iterator */
        __host__ __device__ Iterator begin() noexcept { return Iterator(data()); }
        __host__ __device__ ConstIterator begin() const noexcept { return cbegin(); }
        __host__ __device__ ConstIterator cbegin() const noexcept { return ConstIterator(data()); }
        __host__ __device__ Iterator end() noexcept { return Iterator(data() + Base::getDerived().getLength()); }
        __host__ __device__ ConstIterator end() const noexcept { return cend(); }
        __host__ __device__ ConstIterator cend() const noexcept { return ConstIterator(data() + Base::getDerived().getLength()); }
        __host__ __device__ ReverseIterator rbegin() noexcept { return ReverseIterator(data() + Base::getDerived().getLength() - 1); }
        __host__ __device__ ConstReverseIterator rbegin() const noexcept { return crbegin(); }
        __host__ __device__ ConstReverseIterator crbegin() const noexcept { return ConstReverseIterator(data() + Base::getDerived().getLength() - 1); }
        __host__ __device__ ReverseIterator rend() noexcept { return ReverseIterator(data() - 1); }
        __host__ __device__ ConstReverseIterator rend() const noexcept { return crend(); }
        __host__ __device__ ConstReverseIterator crend() const noexcept { return ConstReverseIterator(data() - 1); }
        /* Getters */
        [[nodiscard]] __host__ __device__ bool empty() const { return Base::getDerived().getLength() == 0; }
        [[nodiscard]] __host__ __device__ PointerType data() noexcept { return Base::getDerived().data(); }
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept { return Base::getDerived().data(); }
    protected:
        __host__ __device__ ArrayBase() = default;
        __host__ __device__ ArrayBase(const ArrayBase&) = default;
        __host__ __device__ ArrayBase(ArrayBase&&) noexcept = default;
        __host__ __device__ ~ArrayBase() = default;
        /* Operators */
        ArrayBase& operator=(const ArrayBase&) = default;
        ArrayBase& operator=(ArrayBase&&) noexcept = default;
    };
}

#include "ArrayBaseImpl.h"
