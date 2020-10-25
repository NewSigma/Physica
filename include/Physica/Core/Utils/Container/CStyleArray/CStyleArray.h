/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_CSTYLEARRAY_H
#define PHYSICA_CSTYLEARRAY_H

#include <cstdlib>
#include <qglobal.h>
#include "AbstractCStyleArray.h"

namespace Physica::Core {
    enum { Dynamic = 0 };
    /*!
     * This class is a wrapper of of array that is allocated by malloc and whose elements is allocated by placement new.
     * This class is designed to avoid passing incorrect arguments to classes that use c-style array only.
     * (Such as \class Matrix and \class Vector.
     * If we pass a array allocated by new to \class Matrix or \class Vector,
     * a memory leak or double free will happen.)
     * Different from containers in other libraries, this class does not call the default constructor of \type T,
     * but directly copy or move existent objects.
     *
     * Note:
     * If \T is a complex type, \T must have its copy and move constructors defined.
     *
     * Unfinished:
     * Copy, move constructors and assign operators maybe able to accept different specializations.
     *
     * Optimize:
     * 1. Use more effective allocate strategy to avoid reallocate.
     * 2. Use the end pointer of arr instead of length may improve performance.
     */
    template<class T, size_t capacity>
    class CStyleArray : public AbstractCStyleArray<T> {
        typedef AbstractCStyleArray<T> Base;
        using Base::arr;
        using Base::length;
    public:
        inline CStyleArray();
        inline explicit CStyleArray(size_t length);
        CStyleArray(std::initializer_list<T> list);
        inline CStyleArray(const CStyleArray& array);
        inline CStyleArray(CStyleArray&& array) noexcept;
        ~CStyleArray() = default;
        /* Operators */
        CStyleArray& operator=(const CStyleArray& array);
        CStyleArray& operator=(CStyleArray&& array) noexcept;
        CStyleArray& operator<<(const T& t) { append(t); return *this; }
        CStyleArray& operator<<(T&& t) { append(std::move(t)); return *this; }
        CStyleArray& operator<<(const CStyleArray& t) { append(t); return *this; }
        CStyleArray& operator<<(CStyleArray&& t) { append(std::move(t)); return *this; }
        /* Helpers */
        CStyleArray<T, Dynamic> subArray(size_t from, size_t to);
        CStyleArray<T, Dynamic> subArray(size_t from) { return subArray(from, Base::getLength()); }
        CStyleArray<T, Dynamic> cut(size_t from);
        void append(const T& t);
        void append(T&& t);
        void append(const CStyleArray& t);
        void append(CStyleArray&& t);
        inline void grow(const T& t);
        inline void grow(T&& t);
        void removeAt(size_t index);
        /*
         * Swap() is not provided, using the swap() of father class instead. That is enough.
         */
        //void swap(CStyleArray& array) noexcept;
        /* Getters */
        [[nodiscard]] constexpr static size_t getCapacity() { return capacity; }
    };

    template<class T>
    class CStyleArray<T, Dynamic> : public AbstractCStyleArray<T> {
        typedef AbstractCStyleArray<T> Base;
        using Base::arr;
        using Base::length;

        size_t capacity;
    public:
        inline CStyleArray();
        inline explicit CStyleArray(size_t length);
        CStyleArray(std::initializer_list<T> list);
        inline CStyleArray(const CStyleArray& array);
        inline CStyleArray(CStyleArray&& array) noexcept;
        ~CStyleArray() = default;
        /* Operators */
        CStyleArray& operator=(const CStyleArray& array);
        CStyleArray& operator=(CStyleArray&& array) noexcept;
        bool operator==(const CStyleArray& array) const;
        bool operator!=(const CStyleArray& array) const { return !(operator==(array)); }
        CStyleArray& operator<<(const T& t) { append(t); return *this; }
        CStyleArray& operator<<(T&& t) { append(std::move(t)); return *this; }
        CStyleArray& operator<<(const CStyleArray& t) { append(t); return *this; }
        CStyleArray& operator<<(CStyleArray&& t) { append(std::move(t)); return *this; }
        /* Helpers */
        CStyleArray<T, Dynamic> subArray(size_t from, size_t to);
        CStyleArray<T, Dynamic> subArray(size_t from) { return subArray(from, length); }
        CStyleArray<T, Dynamic> cut(size_t from);
        void append(const T& t);
        void append(T&& t);
        void append(const CStyleArray& t);
        void append(CStyleArray&& t);
        inline void grow(const T& t);
        inline void grow(T&& t);
        void resize(size_t size);
        void squeeze();
        void increase(size_t size);
        void decrease(size_t size);
        void removeAt(size_t index);
        void swap(CStyleArray& array) noexcept;
        /* Getters */
        [[nodiscard]] size_t getCapacity() const { return capacity; }
    };

    template<class T, size_t capacity>
    inline void swap(CStyleArray<T, capacity>& array1, CStyleArray<T, capacity>& array2) noexcept {
        array1.swap(array2);
    }
}

#include "CStyleArrayImpl.h"

#endif