/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CSTYLEARRAY_H
#define PHYSICA_CSTYLEARRAY_H

#include <cstdlib>
#include <qglobal.h>
#include "AbstractArray.h"

namespace Physica::Core {
    /*!
     * This class is a wrapper of of array that is AbstractCStyleArray<T>::allocated by malloc
     * , whose elements is AbstractCStyleArray<T>::allocated by placement new.
     * This class is designed to avoid passing incorrect arguments to classes that
     * use c-style array only.
     * (Such as \class Matrix and \class Vector.
     * If we pass a array AbstractCStyleArray<T>::allocated by new to \class Matrix or \class Vector,
     * a memory leak or double free will happen.)
     * Different from \class QVector in Qt, this class does not call the default constructor of \type T,
     * but directly copy or move existent elements.
     *
     * Note:
     * If \T is a complex type, \T must have its copy and move constructors defined.
     *
     * Optimize:
     * Copy, move constructors and assign operators maybe able to accept different specializations.
     */
    template<class T, size_t length, size_t capacity = length>
    class CStyleArray : public AbstractArray<T> {
    protected:
        static_assert(length == capacity && length != 0, "Length must not less than capacity.");
        using AbstractArray<T>::arr;
    public:
        inline CStyleArray();
        inline CStyleArray(const CStyleArray<T, length, capacity>& array);
        inline CStyleArray(CStyleArray<T, length, capacity>&& array) noexcept;
        ~CStyleArray();
        /* Operators */
        CStyleArray<T, length, capacity>& operator=(const CStyleArray<T, length, capacity>& array);
        inline CStyleArray<T, length, capacity>& operator=(CStyleArray<T, length, capacity>&& array) noexcept;
        T& operator[](size_t index) { Q_ASSERT(index < length); return arr[index]; }
        const T& operator[](size_t index) const { Q_ASSERT(index < length); return arr[index]; }
        /* Helpers */
        inline void swap(CStyleArray<T, length, capacity>& array) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const { return length; }
        [[nodiscard]] size_t getCapacity() const { return capacity; }
    };

    template<class T, size_t capacity>
    class CStyleArray<T, Dynamic, capacity> : public AbstractArray<T> {
    protected:
        using AbstractArray<T>::arr;
        size_t length;
    public:
        inline CStyleArray();
        inline explicit CStyleArray(size_t length);
        inline CStyleArray(const CStyleArray<T, Dynamic, capacity>& array);
        inline CStyleArray(CStyleArray<T, Dynamic, capacity>&& array) noexcept;
        ~CStyleArray();
        /* Operators */
        CStyleArray<T, Dynamic, capacity>& operator=(const CStyleArray<T, Dynamic, capacity>& array);
        CStyleArray<T, Dynamic, capacity>& operator=(CStyleArray<T, Dynamic, capacity>&& array) noexcept;
        T& operator[](size_t index) { Q_ASSERT(index < length); return arr[index]; }
        const T& operator[](size_t index) const { Q_ASSERT(index < length); return arr[index]; }
        /* Helpers */
        T cutLast();
        void swap(CStyleArray<T, Dynamic, capacity>& array) noexcept;
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
        //!
        void setLength(size_t size) { Q_ASSERT(length <= size); length = size; }
    };

    template<class T>
    class CStyleArray<T, Dynamic, Dynamic> : public AbstractArray<T> {
    protected:
        using AbstractArray<T>::arr;
        size_t length;
        size_t capacity;
    public:
        inline CStyleArray();
        inline explicit CStyleArray(size_t capacity);
        inline CStyleArray(size_t length, size_t capacity);
        inline CStyleArray(const CStyleArray<T, Dynamic, Dynamic>& array);
        inline CStyleArray(CStyleArray<T, Dynamic, Dynamic>&& array) noexcept;
        ~CStyleArray();
        /* Operators */
        CStyleArray<T, Dynamic, Dynamic>& operator=(const CStyleArray<T, Dynamic, Dynamic>& array);
        CStyleArray<T, Dynamic, Dynamic>& operator=(CStyleArray<T, Dynamic, Dynamic>&& array) noexcept;
        T& operator[](size_t index) { Q_ASSERT(index < length); return arr[index]; }
        const T& operator[](size_t index) const { Q_ASSERT(index < length); return arr[index]; }
        CStyleArray<T, Dynamic, Dynamic>& operator<<(const T& t) { append(t); return *this; }
        CStyleArray<T, Dynamic, Dynamic>& operator<<(T&& t) { append(std::move(t)); return *this; }
        /* Helpers */
        CStyleArray<T, Dynamic, Dynamic> subArray(size_t from, size_t to);
        CStyleArray<T, Dynamic, Dynamic> subArray(size_t from) { return subArray(from, length); }
        CStyleArray<T, Dynamic, Dynamic> cut(size_t from);
        void append(const T& t);
        void append(T&& t);
        void append(const CStyleArray<T, Dynamic, Dynamic>& t);
        void append(CStyleArray<T, Dynamic, Dynamic>&& t);
        inline void grow(const T& t);
        inline void grow(T&& t);
        void resize(size_t size);
        void increase(size_t size);
        void decrease(size_t size);
        void squeeze();
        T cutLast();
        void swap(CStyleArray<T, Dynamic, Dynamic>& array) noexcept;
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
        //!
        void setLength(size_t size) { Q_ASSERT(length <= size); length = size; }
    };
}

#include "CStyleArrayImpl.h"

#endif