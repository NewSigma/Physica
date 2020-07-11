/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CSTYLEARRAY_H
#define PHYSICA_CSTYLEARRAY_H

#include <cstdlib>
#include <qglobal.h>

namespace Physica::Core {
    enum { Dynamic = 0 };
    /*!
     * This class is a wrapper of of array that is allocated by malloc and whose elements is allocated by placement new.
     * This class is designed to avoid passing incorrect arguments to classes that use c-style array only.
     * (Such as \class Matrix and \class Vector.
     * If we pass a array CStyleArrayData<T>::allocated by new to \class Matrix or \class Vector,
     * a memory leak or double free will happen.)
     * Different from containers in other libraries, this class does not call the default constructor of \type T,
     * but directly copy or move existent elements.
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
    class CStyleArray {
        T* __restrict arr;
        size_t length;
    public:
        inline CStyleArray();
        inline explicit CStyleArray(size_t length);
        inline CStyleArray(const CStyleArray& array);
        inline CStyleArray(CStyleArray&& array) noexcept;
        ~CStyleArray();
        /* Operators */
        CStyleArray& operator=(const CStyleArray& array);
        CStyleArray& operator=(CStyleArray&& array) noexcept;
        T& operator[](size_t index) { Q_ASSERT(index < length); return arr[index]; }
        const T& operator[](size_t index) const { Q_ASSERT(index < length); return arr[index]; }
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
        inline void allocate(const T& t, size_t index);
        inline void allocate(T&& t, size_t index);
        void append(const T& t);
        void append(T&& t);
        void append(const CStyleArray& t);
        void append(CStyleArray&& t);
        inline void grow(const T& t);
        inline void grow(T&& t);
        void removeAt(size_t index);
        T cutLast();
        void swap(CStyleArray& array) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const { return length; }
        [[nodiscard]] constexpr static size_t getCapacity() { return capacity; }
        [[nodiscard]] bool isEmpty() const { return length == 0; }
        /* Setters */
        /*!
         * Low level api. Designed for performance.
         * \size must larger than current length. Because we can not delete the elements we do not need if not.
         * Elements between old length and \size have not allocated. DO NOT try to visit them.
         */
        void setLength(size_t size) { Q_ASSERT(length <= size); length = size; }
    };

    template<class T>
    class CStyleArray<T, Dynamic> {
        T* __restrict arr;
        size_t length;
        size_t capacity;
    public:
        inline CStyleArray();
        inline explicit CStyleArray(size_t length);
        inline CStyleArray(const CStyleArray& array);
        inline CStyleArray(CStyleArray&& array) noexcept;
        ~CStyleArray();
        /* Operators */
        CStyleArray& operator=(const CStyleArray& array);
        CStyleArray& operator=(CStyleArray&& array) noexcept;
        T& operator[](size_t index) { Q_ASSERT(index < length); return arr[index]; }
        const T& operator[](size_t index) const { Q_ASSERT(index < length); return arr[index]; }
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
        inline void allocate(const T& t, size_t index);
        inline void allocate(T&& t, size_t index);
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
        T cutLast();
        void swap(CStyleArray& array) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const { return length; }
        [[nodiscard]] size_t getCapacity() const { return capacity; }
        [[nodiscard]] bool isEmpty() const { return length == 0; }
        /* Setters */
        /*!
         * Low level api. Designed for performance.
         * \size must larger than current length. Because we can not delete the elements we do not need if not.
         * Elements between old length and \size have not allocated. DO NOT try to visit them.
         */
        void setLength(size_t size) { Q_ASSERT(length <= size); length = size; }
    };

    template<class T, size_t capacity>
    inline void swap(CStyleArray<T, capacity>& array1, CStyleArray<T, capacity>& array2) noexcept {
        array1.swap(array2);
    }
}

#include "CStyleArrayImpl.h"

#endif