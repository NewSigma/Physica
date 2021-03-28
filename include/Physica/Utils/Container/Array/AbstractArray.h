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
#include "Physica/Utils/Template/CRTPBase.h"

namespace Physica::Utils::Internal {
    //Forward declaration
    template<class T> class Traits;
    template<class Derived> class AbstractArray;

    template<class Pointer, class Container>
    class Iterator;
    
    template<class Pointer, class Derived>
    class Iterator<Pointer, AbstractArray<Derived>> {
        Pointer* p;
    public:
        Iterator(const Iterator& ite) : p(ite.p) {}
        ~Iterator() = default;
        /* Operators */
        Iterator& operator=(const Iterator& ite);
        bool operator==(const Iterator& ite) const noexcept { return p == ite.p; }
        bool operator!=(const Iterator& ite) const noexcept { return p != ite.p; }
        Iterator& operator++();
        const Iterator operator++(int);
        Pointer& operator*() { return *p; }
    private:
        explicit Iterator(Pointer* p) : p(p) {}

        friend class AbstractArray<Derived>;
    };

    template<class Pointer, class Container>
    class ReverseIterator;
    
    template<class Pointer, class Derived>
    class ReverseIterator<Pointer, AbstractArray<Derived>> {
        Pointer* p;
    public:
        ReverseIterator(const ReverseIterator& ite) : p(ite.p) {}
        ~ReverseIterator() = default;
        /* Operators */
        ReverseIterator& operator=(const ReverseIterator& ite);
        bool operator==(const ReverseIterator& ite) const noexcept { return p == ite.p; }
        bool operator!=(const ReverseIterator& ite) const noexcept { return p != ite.p; }
        ReverseIterator& operator++();
        const ReverseIterator operator++(int);
        Pointer& operator*() { return *p; }
    private:
        explicit ReverseIterator(Pointer* p) : p(p) {}

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
        using Iterator_ = Iterator<T, AbstractArray<Derived>>;
        using ConstIterator = Iterator<const T, AbstractArray<Derived>>;
        using ReverseIterator_ = ReverseIterator<T, AbstractArray<Derived>>;
        using ConstReverseIterator = ReverseIterator<const T, AbstractArray<Derived>>;

        T* __restrict arr;
    public:
        AbstractArray() = delete;
        /* Operators */
        AbstractArray& operator=(const AbstractArray& array) = delete;
        AbstractArray& operator=(AbstractArray&& array) noexcept = delete;
        T& operator[](size_t index);
        const T& operator[](size_t index) const;
        bool operator==(const AbstractArray& array) const;
        bool operator!=(const AbstractArray& array) const { return !(*this == array); }
        /* Iterator */
        Iterator_ begin() noexcept { return Iterator_(arr); }
        Iterator_ end() noexcept { return Iterator_(arr + getDerived().getLength()); }
        ConstIterator cbegin() const noexcept { return ConstIterator(arr); }
        ConstIterator cend() const noexcept { return ConstIterator(arr + getDerived().getLength()); }
        ReverseIterator_ rbegin() const noexcept { return ReverseIterator_(arr + getDerived().getLength()); }
        ReverseIterator_ rend() const noexcept { return ReverseIterator_(arr - 1); }
        ConstReverseIterator crbegin() const noexcept { return ConstReverseIterator(arr + getDerived().getLength()); }
        ConstReverseIterator crend() const noexcept { return ConstReverseIterator(arr - 1); }
        /* Getters */
        [[nodiscard]] bool empty() const { return getDerived().getLength() == 0; }
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
        /* Getters */
        [[nodiscard]] Derived& getDerived() noexcept { return static_cast<Derived&>(*this); }
        [[nodiscard]] const Derived& getDerived() const noexcept { return static_cast<const Derived&>(*this); }
    };
}

#include "AbstractArrayImpl.h"
