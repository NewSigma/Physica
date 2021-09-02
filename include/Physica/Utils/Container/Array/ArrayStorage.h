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
#include "Physica/Utils/Container/DeivceAllocator.cuh"

namespace Physica::Utils::Internal {
    //Forward declaration
    template<class T> class Traits;
    template<class Derived> class ArrayStorage;

    template<class T, class Container>
    class ContainerIterator;
    
    template<class ValueType, class Derived>
    class ContainerIterator<ValueType, ArrayStorage<Derived>> {
        using PointerType = typename std::add_pointer<ValueType>::type;
        using LValueReferenceType = typename std::add_lvalue_reference<ValueType>::type;

        PointerType p;
    public:
        __host__ __device__ ContainerIterator(const ContainerIterator& ite) : p(ite.p) {}
        ~ContainerIterator() = default;
        /* Operators */
        __host__ __device__ ContainerIterator& operator=(const ContainerIterator& ite);
        __host__ __device__ bool operator==(const ContainerIterator& ite) const noexcept { return p == ite.p; }
        __host__ __device__ bool operator!=(const ContainerIterator& ite) const noexcept { return p != ite.p; }
        __host__ __device__ ContainerIterator& operator++();
        __host__ __device__ const ContainerIterator operator++(int);
        __host__ __device__ LValueReferenceType operator*() const { return *p; }
    private:
        __host__ __device__ explicit ContainerIterator(PointerType p) : p(p) {}

        friend class ArrayStorage<Derived>;
    };

    template<class Pointer, class Container>
    class ReverseContainerIterator;
    
    template<class ValueType, class Derived>
    class ReverseContainerIterator<ValueType, ArrayStorage<Derived>> {
        using PointerType = typename std::add_pointer<ValueType>::type;
        using LValueReferenceType = typename std::add_lvalue_reference<ValueType>::type;

        PointerType p;
    public:
        __host__ __device__ ReverseContainerIterator(const ReverseContainerIterator& ite) : p(ite.p) {}
        ~ReverseContainerIterator() = default;
        /* Operators */
        __host__ __device__ ReverseContainerIterator& operator=(const ReverseContainerIterator& ite);
        __host__ __device__ bool operator==(const ReverseContainerIterator& ite) const noexcept { return p == ite.p; }
        __host__ __device__ bool operator!=(const ReverseContainerIterator& ite) const noexcept { return p != ite.p; }
        __host__ __device__ ReverseContainerIterator& operator++();
        __host__ __device__ const ReverseContainerIterator operator++(int);
        __host__ __device__ LValueReferenceType operator*() const { return *p; }
    private:
        __host__ __device__ explicit ReverseContainerIterator(PointerType p) : p(p) {}

        friend class ArrayStorage<Derived>;
    };
    /**
     * Public parts among specializations of \class Array.
     */
    template<class Derived>
    class ArrayStorage : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using allocator_type = typename Traits<Derived>::AllocatorType;
        using AllocatorTraits = std::allocator_traits<allocator_type>;
        using ValueType = typename AllocatorTraits::value_type;
        using PointerType = typename AllocatorTraits::pointer;
        using LValueReferenceType = typename AllocatorTraits::lvalue_reference;
        using ConstLValueReferenceType = typename AllocatorTraits::const_lvalue_reference;
        using RValueReferenceType = typename AllocatorTraits::rvalue_reference;
    protected:
        PointerType arr;
        allocator_type alloc;
    public:
        using Iterator = ContainerIterator<ValueType, ArrayStorage<Derived>>;
        using ConstIterator = ContainerIterator<const ValueType, ArrayStorage<Derived>>;
        using ReverseIterator = ReverseContainerIterator<ValueType, ArrayStorage<Derived>>;
        using ConstReverseIterator = ReverseContainerIterator<const ValueType, ArrayStorage<Derived>>;
    public:
        ArrayStorage() = delete;
        /* Operators */
        ArrayStorage& operator=(const ArrayStorage& array) = delete;
        ArrayStorage& operator=(ArrayStorage&& array) noexcept = delete;
        [[nodiscard]] __host__ __device__ LValueReferenceType operator[](size_t index);
        [[nodiscard]] __host__ __device__ ConstLValueReferenceType operator[](size_t index) const;
        bool operator==(const ArrayStorage& array) const;
        __host__ __device__ bool operator!=(const ArrayStorage& array) const { return !(*this == array); }
        /* Iterator */
        __host__ __device__ Iterator begin() noexcept { return Iterator(arr); }
        __host__ __device__ Iterator end() noexcept { return Iterator(arr + Base::getDerived().getLength()); }
        __host__ __device__ ConstIterator cbegin() const noexcept { return ConstIterator(arr); }
        __host__ __device__ ConstIterator cend() const noexcept { return ConstIterator(arr + Base::getDerived().getLength()); }
        __host__ __device__ ReverseIterator rbegin() const noexcept { return ReverseIterator(arr + Base::getDerived().getLength()); }
        __host__ __device__ ReverseIterator rend() const noexcept { return ReverseIterator(arr - 1); }
        __host__ __device__ ConstReverseIterator crbegin() const noexcept { return ConstReverseIterator(arr + Base::getDerived().getLength()); }
        __host__ __device__ ConstReverseIterator crend() const noexcept { return ConstReverseIterator(arr - 1); }
        /* Getters */
        [[nodiscard]] __host__ __device__ bool empty() const { return Base::getDerived().getLength() == 0; }
        [[nodiscard]] __host__ __device__ PointerType data() noexcept { return arr; }
        [[nodiscard]] __host__ __device__ const PointerType data() const noexcept { return arr; }
        [[nodiscard]] __host__ __device__ allocator_type get_allocator() const noexcept { return alloc; }
    protected:
        __host__ __device__ explicit ArrayStorage(size_t capacity);
        __host__ __device__ explicit ArrayStorage(PointerType arr_);
        __host__ __device__ ArrayStorage(const ArrayStorage& array);
        __host__ __device__ ArrayStorage(ArrayStorage&& array) noexcept;
        __host__ __device__ ~ArrayStorage();
        /* Helpers */
        __host__ __device__ inline void swap(ArrayStorage& array) noexcept;
    };
}

#include "ArrayStorageImpl.h"
