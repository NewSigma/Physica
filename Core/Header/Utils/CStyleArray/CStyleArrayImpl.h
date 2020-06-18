/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CSTYLEARRAYIMPL_H
#define PHYSICA_CSTYLEARRAYIMPL_H

#ifndef PHYSICA_CSTYLEARRAY_H
    #include "CStyleArray.h"
#endif

namespace Physica::Core {
    ///////////////////////////////////Implementation of CStyleArray<T, length, capacity> //////////////////////////////////
    //!Should be initialized using allocate() before access the elements.
    template<class T, size_t length, size_t capacity>
    inline CStyleArray<T, length, capacity>::CStyleArray() : AbstractArray<T>(capacity) { Q_UNUSED(length) }

    template<class T, size_t length, size_t capacity>
    inline CStyleArray<T, length, capacity>::CStyleArray(const CStyleArray<T, length, capacity>& array)
            : AbstractArray<T>(capacity) {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                new (arr + i) T(array.arr[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T, size_t length, size_t capacity>
    inline CStyleArray<T, length, capacity>::CStyleArray(CStyleArray<T, length, capacity>&& array) noexcept
            : AbstractArray<T>(std::move(array)) {}

    template<class T, size_t length, size_t capacity>
    CStyleArray<T, length, capacity>::~CStyleArray() {
        Q_UNUSED(capacity)
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
        free(arr);
    }

    template<class T, size_t length, size_t capacity>
    CStyleArray<T, length, capacity>&
    CStyleArray<T, length, capacity>::operator=(const CStyleArray<T, length, capacity>& array) {
        if(this != &array) {
            this->~CStyleArray();
            if(QTypeInfo<T>::isComplex)
                for(size_t i = 0; i < length; ++i)
                    new (arr + i) T(array.arr[i]);
            else
                memcpy(arr, array.arr, length * sizeof(T));
        }
        return *this;
    }

    template<class T, size_t length, size_t capacity>
    inline CStyleArray<T, length, capacity>&
    CStyleArray<T, length, capacity>::operator=(CStyleArray<T, length, capacity>&& array) noexcept {
        this->~CStyleArray();
        AbstractArray<T>::operator=(std::move(array));
        return *this;
    }

    template<class T, size_t length, size_t capacity>
    inline void CStyleArray<T, length, capacity>::swap(CStyleArray<T, length, capacity>& array) noexcept {
        AbstractArray<T>::swap(array);
    }

    template<class T, size_t length, size_t capacity>
    inline void swap(CStyleArray<T, length, capacity>& array1, CStyleArray<T, length, capacity>& array2) noexcept {
        array1.swap(array2);
    }
    ///////////////////////////////////Implementation of CStyleArray<T, Dynamic, capacity> //////////////////////////////////
    template<class T, size_t capacity>
    inline CStyleArray<T, Dynamic, capacity>::CStyleArray() : AbstractArray<T>(capacity), length(0) {}

    template<class T, size_t capacity>
    inline CStyleArray<T, Dynamic, capacity>::CStyleArray(size_t length) : AbstractArray<T>(capacity), length(length) {}

    template<class T, size_t capacity>
    inline CStyleArray<T, Dynamic, capacity>::CStyleArray(const CStyleArray<T, Dynamic, capacity>& array)
            : AbstractArray<T>(capacity), length(array.length) {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                new (arr + i) T(array.arr[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T, size_t capacity>
    inline CStyleArray<T, Dynamic, capacity>::CStyleArray(CStyleArray<T, Dynamic, capacity>&& array) noexcept
            : AbstractArray<T>(std::move(array)), length(array.length) {}

    template<class T, size_t capacity>
    CStyleArray<T, Dynamic, capacity>::~CStyleArray() {
        Q_UNUSED(capacity)
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
        free(arr);
    }

    template<class T, size_t capacity>
    CStyleArray<T, Dynamic, capacity>&
    CStyleArray<T, Dynamic, capacity>::operator=(const CStyleArray<T, Dynamic, capacity>& array) {
        if(this != &array) {
            this->~CStyleArray();
            length = array.length;
            if(QTypeInfo<T>::isComplex)
                for(size_t i = 0; i < length; ++i)
                    new (arr + i) T(array.arr[i]);
            else
                memcpy(arr, array.arr, length * sizeof(T));
        }
        return *this;
    }

    template<class T, size_t capacity>
    CStyleArray<T, Dynamic, capacity>&
    CStyleArray<T, Dynamic, capacity>::operator=(CStyleArray<T, Dynamic, capacity>&& array) noexcept {
        this->~CStyleArray();
        AbstractArray<T>::operator=(std::move(array));
        length = array.length;
        array.length = 0;
        return *this;
    }

    template<class T, size_t capacity>
    T CStyleArray<T, Dynamic, capacity>::cutLast() {
        Q_ASSERT(length > 0);
        if(QTypeInfo<T>::isComplex)
            return T(std::move(arr[(length--) - 1]));
        else
            return arr[(length--) - 1];
    }

    template<class T, size_t capacity>
    void CStyleArray<T, Dynamic, capacity>::swap(CStyleArray<T, Dynamic, capacity>& array) noexcept {
        AbstractArray<T>::swap(array);
        std::swap(length, array.length);
    }

    template<class T, size_t capacity>
    inline void swap(CStyleArray<T, Dynamic, capacity>& array1, CStyleArray<T, Dynamic, capacity>& array2) noexcept { array1.swap(array2); }
    ///////////////////////////////////Implementation of CStyleArray<T, Dynamic, Dynamic> ///////////////////////////////
    template<class T>
    inline CStyleArray<T, Dynamic, Dynamic>::CStyleArray() : AbstractArray<T>(), length(0), capacity(0) {}

    template<class T>
    inline CStyleArray<T, Dynamic, Dynamic>::CStyleArray(size_t capacity)
            : AbstractArray<T>(capacity), length(0), capacity(capacity) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not AbstractCStyleArray<T>::allocated. DO NOT try to visit them.
     */
    template<class T>
    inline CStyleArray<T, Dynamic, Dynamic>::CStyleArray(size_t length, size_t capacity)
            : AbstractArray<T>(capacity), length(length), capacity(capacity) {}

    template<class T>
    inline CStyleArray<T, Dynamic, Dynamic>::CStyleArray(const CStyleArray<T, Dynamic, Dynamic>& array)
            : AbstractArray<T>(array.capacity), length(array.length), capacity(array.capacity) {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                new (arr + i) T(array.arr[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T>
    inline CStyleArray<T, Dynamic, Dynamic>::CStyleArray(CStyleArray<T, Dynamic, Dynamic>&& array) noexcept
            : AbstractArray<T>(std::move(array)), length(array.length), capacity(array.capacity) {}

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::~CStyleArray() {
        Q_UNUSED(capacity)
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
        free(arr);
    }

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>& CStyleArray<T, Dynamic, Dynamic>::operator=(const CStyleArray<T, Dynamic, Dynamic>& array) {
        if(this != &array) {
            this->~CStyleArray();
            length = array.length;
            capacity = array.capacity;
            if(QTypeInfo<T>::isComplex)
                for(size_t i = 0; i < length; ++i)
                    new (arr + i) T(array.arr[i]);
            else
                memcpy(arr, array.arr, length * sizeof(T));
        }
        return *this;
    }

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>&
    CStyleArray<T, Dynamic, Dynamic>::operator=(CStyleArray<T, Dynamic, Dynamic>&& array) noexcept {
        this->~CStyleArray();
        length = array.length;
        array.length = 0;
        capacity = array.capacity;
        AbstractArray<T>::operator=(std::move(array));
        return *this;
    }
    /*!
     * Return the sub array of current array. From is included and to is excluded.
     */
    template<class T>
    CStyleArray<T, Dynamic, Dynamic> CStyleArray<T, Dynamic, Dynamic>::subArray(size_t from, size_t to) {
        Q_ASSERT(from < to && to <= length);
        const auto result_length = to - from;
        CStyleArray<T, Dynamic, Dynamic> result(result_length);
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < result_length; ++i, ++from)
                new (result.arr + i) T(arr[from]);
        else
            memcpy(result.arr, arr, result_length * sizeof(T));
    }
    //!\from is included
    template<class T>
    CStyleArray<T, Dynamic, Dynamic> CStyleArray<T, Dynamic, Dynamic>::cut(size_t from) {
        Q_ASSERT(from < length);
        auto result_length = length - from;
        CStyleArray<T, Dynamic, Dynamic> result(result_length);
        length = from;
        for(size_t i = 0; from < length; ++from, ++i)
            result.grow(std::move(arr[from]));
        return result;
    }
    /*!
     * Wrap structure: append() <- grow() <- AbstractCStyleArray<T>::allocate()
     */
    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(const T& t) {
        if(length == capacity)
            increase(capacity + 1);
        grow(t);
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(T&& t) {
        if(length == capacity)
            increase(capacity + 1);
        grow(std::move(t));
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(const CStyleArray<T, Dynamic, Dynamic>& t) {
        const auto new_length = length + t.length;
        if(new_length > capacity)
            resize(new_length);
        for(size_t i = 0; i < t.length; ++i, ++length)
            AbstractArray<T>::allocate(t[i], length);
        length = new_length;
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(CStyleArray<T, Dynamic, Dynamic>&& t) {
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
    inline void CStyleArray<T, Dynamic, Dynamic>::grow(const T& t) {
        Q_ASSERT(length < capacity);
        AbstractArray<T>::allocate(t, length++);
    }

    template<class T>
    inline void CStyleArray<T, Dynamic, Dynamic>::grow(T&& t) {
        Q_ASSERT(length < capacity);
        AbstractArray<T>::allocate(std::move(t), length++);
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::resize(size_t size) {
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
    void CStyleArray<T, Dynamic, Dynamic>::increase(size_t size) {
        Q_ASSERT(size >= capacity);
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        capacity = size;
    }
    /*!
     * Decrease the capacity.
     * This function can be used when you are sure the new \size is shorter than the old capacity.
     */
    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::decrease(size_t size) {
        Q_ASSERT(size <= capacity);
        if(QTypeInfo<T>::isComplex) {
            for(size_t i = size; size < length; ++i)
                (arr + i)->~T();
        }
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        length = capacity = size;
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::squeeze() {
        arr = reinterpret_cast<T*>(realloc(arr, length * sizeof(T)));
        capacity = length;
    }

    template<class T>
    T CStyleArray<T, Dynamic, Dynamic>::cutLast() {
        Q_ASSERT(length > 0);
        if(QTypeInfo<T>::isComplex)
            return T(std::move(arr[(length--) - 1]));
        else
            return arr[(length--) - 1];
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::swap(CStyleArray<T, Dynamic, Dynamic>& array) noexcept {
        std::swap(length, array.length);
        std::swap(capacity, array.capacity);
        AbstractArray<T>::swap(array);
    }

    template<class T>
    inline void swap(CStyleArray<T, Dynamic, Dynamic>& array1, CStyleArray<T, Dynamic, Dynamic>& array2) noexcept {
        array1.swap(array2);
    }
}

#endif