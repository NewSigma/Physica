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

#include <cassert>
#include <cstdlib>
#include <qglobal.h>
#include "AbstractCStyleArrayWithLength.h"

namespace Physica::Utils {
    constexpr size_t Dynamic = 0;
    /*!
     * This class is a wrapper of of array that is allocated by malloc and whose elements is allocated by placement new.
     * This class is designed to avoid passing incorrect arguments to classes that use c-style array only.
     * (Such as \class Matrix and \class Vector.
     * If we pass a array allocated by new to \class Matrix or \class Vector,
     * a memory leak or double free will happen.)
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

    namespace Intenal {
        template<class T, size_t Length, size_t Capacity>
        class Traits<CStyleArray<T, Length, Capacity>> {
        public:
            using ElementType = T;
            constexpr static size_t ArrayLength = Length;
            constexpr static size_t ArrayCapacity = Capacity;
        };
    }

    template<class T, size_t Length, size_t Capacity>
    class CStyleArray : public Intenal::AbstractCStyleArray<CStyleArray<T, Length, Capacity>> {
        static_assert(Length == Capacity, "Capacity of fixed array must equals to Length.");
    private:
        using Base = Intenal::AbstractCStyleArray<CStyleArray<T, Length, Capacity>>;
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
        CStyleArray<T, Dynamic, Dynamic> subArray(size_t from, size_t to);
        CStyleArray<T, Dynamic, Dynamic> subArray(size_t from) { return subArray(from, Length); }
        CStyleArray<T, Dynamic, Dynamic> cut(size_t from);
        void allocate(const T& t, size_t index) { assert(false); } //Never call it, for the convience of implement templates.
        void allocate(T&& t, size_t index) { assert(false); } //Never call it, for the convience of implement templates.
        void swap(CStyleArray& array) noexcept { Base::swap(array); }
        /* Getters */
        [[nodiscard]] constexpr static size_t getLength() { return Length; }
        [[nodiscard]] constexpr static size_t getCapacity() { return Capacity; }
    };

    template<class T, size_t Capacity>
    class CStyleArray<T, Dynamic, Capacity>
        : public Intenal::AbstractCStyleArrayWithLength<CStyleArray<T, Dynamic, Capacity>> {
    private:
        using Base = Intenal::AbstractCStyleArrayWithLength<CStyleArray<T, Dynamic, Capacity>>;
        using Base::length;
        using Base::arr;
    public:
        CStyleArray();
        explicit CStyleArray(size_t length);
        CStyleArray(std::initializer_list<T> list);
        CStyleArray(const CStyleArray& array);
        CStyleArray(CStyleArray&& array) noexcept;
        ~CStyleArray() = default;
        /* Operators */
        CStyleArray& operator=(const CStyleArray& array);
        CStyleArray& operator=(CStyleArray&& array) noexcept;
        /* Helpers */
        CStyleArray<T, Dynamic, Dynamic> subArray(size_t from, size_t to);
        CStyleArray<T, Dynamic, Dynamic> subArray(size_t from) { return subArray(from, length); }
        CStyleArray<T, Dynamic, Dynamic> cut(size_t from);
        inline void append(const T& t);
        inline void append(T&& t);
        void append(const CStyleArray& t);
        void append(CStyleArray&& t);
        void swap(CStyleArray& array) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return length; }
        [[nodiscard]] constexpr static size_t getCapacity() { return Capacity; }
    };

    template<class T>
    class CStyleArray<T, Dynamic, Dynamic>
        : public Intenal::AbstractCStyleArrayWithLength<CStyleArray<T, Dynamic, Dynamic>> {
    private:
        using Base = Intenal::AbstractCStyleArrayWithLength<CStyleArray<T, Dynamic, Dynamic>>;
        using Base::length;
        using Base::arr;
    protected:
        size_t capacity;
    public:
        CStyleArray();
        explicit CStyleArray(size_t length_);
        CStyleArray(size_t length_, size_t capacity_);
        CStyleArray(std::initializer_list<T> list);
        CStyleArray(const CStyleArray& array);
        CStyleArray(CStyleArray&& array) noexcept;
        ~CStyleArray() = default;
        /* Operators */
        CStyleArray& operator=(const CStyleArray& array);
        CStyleArray& operator=(CStyleArray&& array) noexcept;
        /* Helpers */
        void append(const T& t);
        void append(T&& t);
        void append(const CStyleArray& t);
        void append(CStyleArray&& t);
        void resize(size_t size);
        void squeeze();
        void increase(size_t size);
        void decrease(size_t size);
        void swap(CStyleArray& array) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return length; }
        [[nodiscard]] size_t getCapacity() const noexcept { return capacity; }
    };

    template<class T, size_t Length, size_t Capacity>
    inline void swap(CStyleArray<T, Length, Capacity>& array1, CStyleArray<T, Length, Capacity>& array2) noexcept {
        array1.swap(array2);
    }
}

#include "CStyleArrayImpl.h"
