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
    template<class Derived> class ArrayStorage;

    template<class T, class Container>
    class ContainerIterator;
    
    template<class ValueType, class Derived>
    class ContainerIterator<ValueType, ArrayStorage<Derived>> {
        using PointerType = typename std::add_pointer<ValueType>::type;
        using LValueReferenceType = typename std::add_lvalue_reference<ValueType>::type;

        PointerType p;
    public:
        ContainerIterator(const ContainerIterator& ite) : p(ite.p) {}
        ~ContainerIterator() = default;
        /* Operators */
        ContainerIterator& operator=(const ContainerIterator& ite);
        bool operator==(const ContainerIterator& ite) const noexcept { return p == ite.p; }
        bool operator!=(const ContainerIterator& ite) const noexcept { return p != ite.p; }
        ContainerIterator& operator++();
        const ContainerIterator operator++(int);
        LValueReferenceType operator*() const { return *p; }
    private:
        explicit ContainerIterator(PointerType p) : p(p) {}

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
        ReverseContainerIterator(const ReverseContainerIterator& ite) : p(ite.p) {}
        ~ReverseContainerIterator() = default;
        /* Operators */
        ReverseContainerIterator& operator=(const ReverseContainerIterator& ite);
        bool operator==(const ReverseContainerIterator& ite) const noexcept { return p == ite.p; }
        bool operator!=(const ReverseContainerIterator& ite) const noexcept { return p != ite.p; }
        ReverseContainerIterator& operator++();
        const ReverseContainerIterator operator++(int);
        LValueReferenceType operator*() const { return *p; }
    private:
        explicit ReverseContainerIterator(PointerType p) : p(p) {}

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
        [[nodiscard]] LValueReferenceType operator[](size_t index);
        [[nodiscard]] ConstLValueReferenceType operator[](size_t index) const;
        bool operator==(const ArrayStorage& array) const;
        bool operator!=(const ArrayStorage& array) const { return !(*this == array); }
        /* Iterator */
        Iterator begin() noexcept { return Iterator(arr); }
        Iterator end() noexcept { return Iterator(arr + Base::getDerived().getLength()); }
        ConstIterator cbegin() const noexcept { return ConstIterator(arr); }
        ConstIterator cend() const noexcept { return ConstIterator(arr + Base::getDerived().getLength()); }
        ReverseIterator rbegin() const noexcept { return ReverseIterator(arr + Base::getDerived().getLength()); }
        ReverseIterator rend() const noexcept { return ReverseIterator(arr - 1); }
        ConstReverseIterator crbegin() const noexcept { return ConstReverseIterator(arr + Base::getDerived().getLength()); }
        ConstReverseIterator crend() const noexcept { return ConstReverseIterator(arr - 1); }
        /* Getters */
        [[nodiscard]] bool empty() const { return Base::getDerived().getLength() == 0; }
        [[nodiscard]] PointerType data() noexcept { return arr; }
        [[nodiscard]] const PointerType data() const noexcept { return arr; }
        [[nodiscard]] allocator_type get_allocator() const noexcept { return alloc; }
    protected:
        explicit ArrayStorage(size_t capacity);
        explicit ArrayStorage(PointerType arr_);
        ArrayStorage(const ArrayStorage& array);
        ArrayStorage(ArrayStorage&& array) noexcept;
        ~ArrayStorage();
        /* Helpers */
        inline void swap(ArrayStorage& array) noexcept;
    };
}

#include "ArrayStorageImpl.h"
