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
#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Utils::Internal {
    //Forward declaration
    template<class T> class Traits;
    template<class Derived> class AbstractArray;

    template<class Pointer, class Container>
    class ContainerIterator;
    
    template<class Pointer, class Derived>
    class ContainerIterator<Pointer, AbstractArray<Derived>> {
        Pointer* p;
    public:
        ContainerIterator(const ContainerIterator& ite) : p(ite.p) {}
        ~ContainerIterator() = default;
        /* Operators */
        ContainerIterator& operator=(const ContainerIterator& ite);
        bool operator==(const ContainerIterator& ite) const noexcept { return p == ite.p; }
        bool operator!=(const ContainerIterator& ite) const noexcept { return p != ite.p; }
        ContainerIterator& operator++();
        const ContainerIterator operator++(int);
        Pointer& operator*() const { return *p; }
    private:
        explicit ContainerIterator(Pointer* p) : p(p) {}

        friend class AbstractArray<Derived>;
    };

    template<class Pointer, class Container>
    class ReverseContainerIterator;
    
    template<class Pointer, class Derived>
    class ReverseContainerIterator<Pointer, AbstractArray<Derived>> {
        Pointer* p;
    public:
        ReverseContainerIterator(const ReverseContainerIterator& ite) : p(ite.p) {}
        ~ReverseContainerIterator() = default;
        /* Operators */
        ReverseContainerIterator& operator=(const ReverseContainerIterator& ite);
        bool operator==(const ReverseContainerIterator& ite) const noexcept { return p == ite.p; }
        bool operator!=(const ReverseContainerIterator& ite) const noexcept { return p != ite.p; }
        ReverseContainerIterator& operator++();
        const ReverseContainerIterator operator++(int);
        Pointer& operator*() const { return *p; }
    private:
        explicit ReverseContainerIterator(Pointer* p) : p(p) {}

        friend class AbstractArray<Derived>;
    };
    /**
     * Public parts among specializations of \class Array.
     */
    template<class Derived>
    class AbstractArray : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    protected:
        using T = typename Traits<Derived>::ElementType;

        T* __restrict arr;
    public:
        using Iterator = ContainerIterator<T, AbstractArray<Derived>>;
        using ConstIterator = ContainerIterator<const T, AbstractArray<Derived>>;
        using ReverseIterator = ReverseContainerIterator<T, AbstractArray<Derived>>;
        using ConstReverseIterator = ReverseContainerIterator<const T, AbstractArray<Derived>>;
    public:
        AbstractArray() = delete;
        /* Operators */
        AbstractArray& operator=(const AbstractArray& array) = delete;
        AbstractArray& operator=(AbstractArray&& array) noexcept = delete;
        [[nodiscard]] T& operator[](size_t index);
        [[nodiscard]] const T& operator[](size_t index) const;
        bool operator==(const AbstractArray& array) const;
        bool operator!=(const AbstractArray& array) const { return !(*this == array); }
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
        [[nodiscard]] T* data() noexcept { return arr; }
        [[nodiscard]] const T* data() const noexcept { return arr; }
    protected:
        explicit AbstractArray(size_t capacity);
        explicit AbstractArray(T* __restrict arr_);
        AbstractArray(const AbstractArray& array);
        AbstractArray(AbstractArray&& array) noexcept;
        ~AbstractArray();
        /* Helpers */
        inline void allocate(const T& t, size_t index);
        inline void allocate(T&& t, size_t index);
        inline void swap(AbstractArray& array) noexcept;
    };
}

#include "AbstractArrayImpl.h"
