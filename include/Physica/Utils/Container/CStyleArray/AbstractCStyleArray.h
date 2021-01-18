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

namespace Physica::Utils::Intenal {
    //Forward declaration
    template<class T> class Traits;
    template<class Derived> class AbstractCStyleArray;

    template<class Pointer, class Container>
    class Iterator;
    
    template<class Pointer, class Derived>
    class Iterator<Pointer, AbstractCStyleArray<Derived>> {
        Pointer* p;
    public:
        Iterator(const Iterator& ite) : p(ite.p) {}
        Iterator(Iterator&& ite) noexcept : p(ite.p) {}
        ~Iterator() = default;
        /* Operators */
        Iterator& operator=(const Iterator& ite);
        Iterator& operator=(Iterator&& ite) noexcept;
        bool operator==(const Iterator& ite) const noexcept { return p == ite.p; }
        bool operator!=(const Iterator& ite) const noexcept { return p != ite.p; }
        Iterator& operator++();
        const Iterator operator++(int);
        Pointer& operator*() { return *p; }
    private:
        explicit Iterator(Pointer* p) : p(p) {}

        friend class AbstractCStyleArray<Derived>;
    };
    /**
     * Public parts among specializations of CStyleArray.
     */
    template<class Derived>
    class AbstractCStyleArray : public Utils::CRTPBase<Derived> {
    protected:
        using T = typename Traits<Derived>::ElementType;
        using Iterator_ = Iterator<T, AbstractCStyleArray<Derived>>;
        using ConstIterator = Iterator<const T, AbstractCStyleArray<Derived>>;
    public:
        T* __restrict arr;
    public:
        AbstractCStyleArray() = delete;
        /* Operators */
        inline T& operator[](size_t index);
        inline const T& operator[](size_t index) const;
        bool operator==(const AbstractCStyleArray& array) const;
        bool operator!=(const AbstractCStyleArray& array) const { return !(operator==(array)); }
        /* Iterator */
        Iterator_ begin() noexcept { return Iterator_(arr); }
        Iterator_ end() noexcept { return Iterator_(arr + getDerived().getLength()); }
        ConstIterator cbegin() const noexcept { return ConstIterator(arr); }
        ConstIterator cend() const noexcept { return ConstIterator(arr + getDerived().getLength()); }
        /* Getters */
        [[nodiscard]] bool empty() const { return getDerived().getLength() == 0; }
    protected:
        explicit AbstractCStyleArray(size_t capacity);
        AbstractCStyleArray(const AbstractCStyleArray& array);
        AbstractCStyleArray(AbstractCStyleArray&& array) noexcept;
        ~AbstractCStyleArray();
        /* Operators */
        AbstractCStyleArray& operator=(const AbstractCStyleArray& array);
        AbstractCStyleArray& operator=(AbstractCStyleArray&& array) noexcept;
        /* Helpers */
        inline void allocate(const T& t, size_t index);
        inline void allocate(T&& t, size_t index);
        inline void swap(AbstractCStyleArray& array) noexcept;
        /* Getters */
        [[nodiscard]] Derived& getDerived() noexcept { return static_cast<Derived&>(*this); }
        [[nodiscard]] const Derived& getDerived() const noexcept { return static_cast<const Derived&>(*this); }
    };
}

#include "AbstractCStyleArrayImpl.h"
