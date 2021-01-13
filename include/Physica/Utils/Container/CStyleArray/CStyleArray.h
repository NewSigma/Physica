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

#include <cstdlib>
#include <qglobal.h>
#include "AbstractCStyleArray.h"

namespace Physica::Utils {
    constexpr size_t Dynamic = 0;
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
    template<class T, size_t Length = Dynamic, size_t Capacity = Length>
    class CStyleArray;

    template<class T, size_t Length, size_t Capacity>
    class CStyleArray<T, Length, Capacity> : public Intenal::AbstractCStyleArray<CStyleArray<T, Length, Capacity>> {
        static_assert(Length == Capacity, "Capacity of fixed array must equals to Length.");
    private:
        using Base = Intenal::AbstractCStyleArray<CStyleArray<T, Length, Capacity>>;
        using Base::Iterator;
        using Base::arr;
    public:
        CStyleArray();
        CStyleArray(std::initializer_list<T> list);
        CStyleArray(const CStyleArray& array);
        CStyleArray(CStyleArray&& array) noexcept;
        ~CStyleArray();
        /* Operators */
        CStyleArray& operator=(const CStyleArray& array);
        CStyleArray& operator=(CStyleArray&& array) noexcept;
        /* Helpers */
        CStyleArray<T, Dynamic> subArray(size_t from, size_t to);
        CStyleArray<T, Dynamic> subArray(size_t from) { return subArray(from, Length); }
        CStyleArray<T, Dynamic> cut(size_t from);
        void swap(CStyleArray& array) noexcept { Base::swap(array); }
        /* Getters */
        [[nodiscard]] constexpr static size_t getLength() { return Length; }
        [[nodiscard]] constexpr static size_t getCapacity() { return Capacity; }
    };

    template<class T, size_t Capacity>
    class CStyleArray<T, Dynamic, Capacity> : public Intenal::AbstractCStyleArray<CStyleArray<T, Dynamic, Capacity>> {
    private:
        using Base = Intenal::AbstractCStyleArray<CStyleArray<T, Dynamic, Capacity>>;
        using Base::Iterator;
        using Base::arr;

        size_t length;
    public:
        CStyleArray();
        explicit CStyleArray(size_t length);
        CStyleArray(std::initializer_list<T> list);
        CStyleArray(const CStyleArray& array);
        CStyleArray(CStyleArray&& array) noexcept;
        ~CStyleArray();
        /* Operators */
        CStyleArray& operator=(const CStyleArray& array);
        CStyleArray& operator=(CStyleArray&& array) noexcept;
        bool operator==(const CStyleArray& array) const;
        bool operator!=(const CStyleArray& array) const { return !(operator==(array)); }
        CStyleArray& operator<<(const T& t) { append(t); return *this; }
        CStyleArray& operator<<(T&& t) { append(std::move(t)); return *this; }
        CStyleArray& operator<<(const CStyleArray& array) { append(array); return *this; }
        CStyleArray& operator<<(CStyleArray&& array) { append(std::move(array)); return *this; }
        /* Helpers */
        T cutLast();
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
        [[nodiscard]] size_t getLength() noexcept const { return length; }
        [[nodiscard]] constexpr static size_t getCapacity() { return Capacity; }
    protected:
        /* Setters */
        /**
         * Low level api. Designed for performance.
         * \size must larger than current length. Because we can not delete the elements we do not need if not.
         * Elements between old length and \size have not allocated. DO NOT try to visit them.
         */
        void setLength(size_t size) { Q_ASSERT(length <= size); length = size; }
    };

    template<class T>
    class CStyleArray<T, Dynamic, Dynamic> : public Intenal::AbstractCStyleArray<CStyleArray<T, Dynamic, Dynamic>> {
    private:
        using Base = Intenal::AbstractCStyleArray<CStyleArray<T, Dynamic, Dynamic>>;
        using Base::Iterator;
        using Base::arr;

        size_t length;
        size_t capacity;
    public:
        CStyleArray();
        explicit CStyleArray(size_t length);
        CStyleArray(std::initializer_list<T> list);
        CStyleArray(const CStyleArray& array);
        CStyleArray(CStyleArray&& array);
        ~CStyleArray();
        /* Operators */
        CStyleArray& operator=(const CStyleArray& array);
        CStyleArray& operator=(CStyleArray&& array) noexcept;
        CStyleArray& operator<<(const T& t) { append(t); return *this; }
        CStyleArray& operator<<(T&& t) { append(std::move(t)); return *this; }
        CStyleArray& operator<<(const CStyleArray& array) { append(array); return *this; }
        CStyleArray& operator<<(CStyleArray&& array) { append(std::move(array)); return *this; }
        /* Helpers */
        T cutLast();
        void append(const T& t);
        void append(T&& t);
        void append(const CStyleArray& t);
        void append(CStyleArray&& t);
        inline void grow(const T& t);
        inline void grow(T&& t);
        void removeAt(size_t index);
        void swap(CStyleArray& array);
        /* Getters */
        [[nodiscard]] size_t getLength() noexcept const { return length; }
        [[nodiscard]] size_t getCapacity() noexcept const { return capacity; }
    protected:
        /* Setters */
        /**
         * Low level api. Designed for performance.
         * \size must larger than current length. Because we can not delete the elements we do not need if not.
         * Elements between old length and \size have not allocated. DO NOT try to visit them.
         */
        void setLength(size_t size) { Q_ASSERT(length <= size); length = size; }
    }

    template<class T, size_t Length, size_t Capacity>
    inline void swap(CStyleArray<T, Length, Capacity>& array1, CStyleArray<T, Length, Capacity>& array2) noexcept {
        array1.swap(array2);
    }
}

#include "CStyleArrayImpl.h"
