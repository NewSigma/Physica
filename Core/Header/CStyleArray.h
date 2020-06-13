/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CSTYLEARRAY_H
#define PHYSICA_CSTYLEARRAY_H

#include <cstdlib>
#include <qglobal.h>

namespace Physica::Core {
    /*!
    * This class is a wrapper of of array that is allocated by malloc, whose elements is allocated
    * by placement new. This class is designed to avoid passing incorrect arguments to classes that
    * use c-style array only.
    * (Such as \class Matrix and \class Vector.
    * If we pass a array allocated by new to \class Matrix or \class Vector,
    * a memory leak or double free will happen.)
    *
    * Note:
    * If \T is a complex type, T must have its copy and move constructors defined.
    *
    * Different from \class QVector in Qt, this class does not call the default constructor of \type T.
    */
    template<class T>
    class CStyleArray {
        T* __restrict arr;
        size_t length;
        size_t capacity;
    public:
        CStyleArray();
        explicit CStyleArray(size_t capacity);
        CStyleArray(size_t length, size_t capacity);
        CStyleArray(const CStyleArray<T>& array);
        CStyleArray(CStyleArray<T>&& array) noexcept;
        ~CStyleArray();
        /* Operators */
        CStyleArray& operator=(const CStyleArray<T>& array);
        CStyleArray& operator=(CStyleArray<T>&& array) noexcept;
        T& operator[](size_t index) { Q_ASSERT(index < length); return arr[index]; }
        const T& operator[](size_t index) const { Q_ASSERT(index < length); return arr[index]; }
        CStyleArray& operator<<(const T& t) { append(t); return *this; }
        CStyleArray& operator<<(T&& t) { append(t); return *this; }
        /* Helpers */
        CStyleArray<T> subArray(size_t from, size_t to);
        CStyleArray<T> subArray(size_t from) { return subArray(from, length); }
        void append(const T& t);
        void append(T&& t);
        void append(const CStyleArray<T>& t);
        void append(CStyleArray<T>&& t);
        inline void grow(const T& t);
        inline void grow(T&& t);
        inline void allocate(const T& t, size_t index);
        inline void allocate(T&& t, size_t index);
        CStyleArray<T> cut(size_t from);
        T cutLast();
        void resize(size_t size);
        void increase(size_t size);
        void decrease(size_t size);
        void squeeze();
        void swap(CStyleArray<T>& array) noexcept;
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
    CStyleArray<T>::CStyleArray() : arr(nullptr), length(0), capacity(0) {}

    template<class T>
    CStyleArray<T>::CStyleArray(size_t capacity)
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))), length(0), capacity(capacity) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not allocated. DO NOT try to visit them.
     */
    template<class T>
    CStyleArray<T>::CStyleArray(size_t length, size_t capacity)
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))), length(length), capacity(capacity) {}

    template<class T>
    CStyleArray<T>::CStyleArray(const CStyleArray<T>& array)
            : arr(reinterpret_cast<T*>(malloc(array.capacity * sizeof(T))))
            , length(array.length), capacity(array.capacity) {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                new (arr + i) T(array.arr[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T>
    CStyleArray<T>::CStyleArray(CStyleArray<T>&& array) noexcept
            : arr(array.arr), length(array.length), capacity(array.capacity) {
        array.arr = nullptr;
        array.length = 0;
    }

    template<class T>
    CStyleArray<T>::~CStyleArray() {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
        free(arr);
    }

    template<class T>
    CStyleArray<T>& CStyleArray<T>::operator=(const CStyleArray<T>& array) {
        if(this == &array)
            return *this;
        this->~CStyleArray();
        length = array.length;
        capacity = array.capacity;
        arr = reinterpret_cast<T*>(malloc(capacity * sizeof(T)));
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                new (arr + i) T(array.arr[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T>
    CStyleArray<T>& CStyleArray<T>::operator=(CStyleArray<T>&& array) noexcept {
        this->~CStyleArray();
        arr = array.arr;
        length = array.length;
        capacity = array.capacity;
        array.arr = nullptr;
        array.length = 0;
    }
    /*!
     * Return the sub array of current array. From is included and to is excluded.
     */
    template<class T>
    CStyleArray<T> CStyleArray<T>::subArray(size_t from, size_t to) {
        Q_ASSERT(from < to && to <= length);
        const auto result_length = to - from;
        CStyleArray<T> result(result_length);
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < result_length; ++i, ++from)
                new (result.arr + i) T(arr[from]);
        else
            memcpy(result.arr, arr, result_length * sizeof(T));
    }
    /*!
     * Wrap structure: append() <- grow() <- allocate()
     */
    template<class T>
    void CStyleArray<T>::append(const T& t) {
        if(length == capacity)
            increase(capacity + 1);
        grow(t);
    }

    template<class T>
    void CStyleArray<T>::append(T&& t) {
        if(length == capacity)
            increase(capacity + 1);
        grow(t);
    }

    template<class T>
    void CStyleArray<T>::append(const CStyleArray<T>& t) {
        const auto new_length = length + t.length;
        if(new_length > capacity)
            resize(new_length);
        for(size_t i = 0; i < t.length; ++i, ++length)
            allocate(t[i], length);
        length = new_length;
    }

    template<class T>
    void CStyleArray<T>::append(CStyleArray<T>&& t) {
        const auto new_length = length + t.length;
        if(new_length > capacity)
            resize(new_length);
        memcpy(arr + length, t.arr, t.length * sizeof(T));
        length = new_length;
        t.arr = nullptr;
        t.length = 0;
    }
    /*!
     * Low level api. Designed for performance.
     * Increase the capacity.
     * This function can be used when you are sure the current capacity is enough.
     */
    template<class T>
    inline void CStyleArray<T>::grow(const T& t) {
        Q_ASSERT(length < capacity);
        allocate(t, length++);
    }

    template<class T>
    inline void CStyleArray<T>::grow(T&& t) {
        Q_ASSERT(length < capacity);
        allocate(std::move(t), length++);
    }
    /*!
     * Low level api. Designed for performance.
     * Simply allocate a \T at position \index.
     * You must ensure position \index is not used or a memory leak will occur.
     */
    template<class T>
    inline void CStyleArray<T>::allocate(const T& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(t);
        else
            *(arr + index) = t;
    }

    template<class T>
    inline void CStyleArray<T>::allocate(T&& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(t);
        else
            *(arr + index) = t;
    }
    //!\from is included
    template<class T>
    CStyleArray<T> CStyleArray<T>::cut(size_t from) {
        Q_ASSERT(from < length);
        auto result_length = length - from;
        CStyleArray<T> result(result_length);
        length = from;
        for(size_t i = 0; from < length; ++from, ++i)
            result.grow(std::move(arr[from]));
        return result;
    }

    template<class T>
    T CStyleArray<T>::cutLast() {
        Q_ASSERT(length > 0);
        if(QTypeInfo<T>::isComplex)
            return T(std::move(arr[(length--) - 1]));
        else
            return arr[(length--) - 1];
    }

    template<class T>
    void CStyleArray<T>::resize(size_t size) {
        if(QTypeInfo<T>::isComplex) {
            if(length > size) {
                for(size_t i = size; size < length; ++i)
                    (arr + i)->~T();
                length = size;
            }
        }
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        capacity = size;
    }
    /*!
     * Increase the capacity.
     * This function can be used when you are sure the new \size is larger than the old capacity.
     */
    template<class T>
    void CStyleArray<T>::increase(size_t size) {
        Q_ASSERT(size >= capacity);
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        capacity = size;
    }
    /*!
     * Decrease the capacity.
     * This function can be used when you are sure the new \size is shorter than the old capacity.
     */
    template<class T>
    void CStyleArray<T>::decrease(size_t size) {
        Q_ASSERT(size <= capacity);
        if(QTypeInfo<T>::isComplex) {
            for(size_t i = size; size < length; ++i)
                (arr + i)->~T();
        }
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        length = capacity = size;
    }

    template<class T>
    void CStyleArray<T>::squeeze() {
        arr = reinterpret_cast<T*>(realloc(arr, length * sizeof(T)));
        capacity = length;
    }

    template<class T>
    void CStyleArray<T>::swap(CStyleArray<T>& array) noexcept {
        auto temp = arr;
        arr = array.arr;
        array.arr = temp;
        auto temp_size = length;
        length = array.length;
        array.length = temp_size;
        temp_size = capacity;
        capacity = array.capacity;
        array.capacity = temp_size;
    }

    template<class T>
    void swap(CStyleArray<T> array1, CStyleArray<T> array2) noexcept { array1.swap(array2); }
}

#endif
