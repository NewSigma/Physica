/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CSTYLEARRAYIMPL_H
#define PHYSICA_CSTYLEARRAYIMPL_H

namespace Physica::Core {
    //////////////////////////////////////////CStyleArray<T, capacity>//////////////////////////////////////////
    template<class T, size_t capacity>
    inline CStyleArray<T, capacity>::CStyleArray()
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))), length(0) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not AbstractCStyleArray<T>::allocated. DO NOT try to visit them.
     */
    template<class T, size_t capacity>
    inline CStyleArray<T, capacity>::CStyleArray(size_t length)
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))), length(length) {}

    template<class T, size_t capacity>
    inline CStyleArray<T, capacity>::CStyleArray(const CStyleArray<T, capacity>& array)
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))), length(array.length) {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                new (arr + i) T(array.arr[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T, size_t capacity>
    inline CStyleArray<T, capacity>::CStyleArray(CStyleArray<T, capacity>&& array) noexcept
            : arr(array.arr), length(array.length) {
        array.arr = nullptr;
        array.length = 0;
    }

    template<class T, size_t capacity>
    CStyleArray<T, capacity>::~CStyleArray() {
        Q_UNUSED(capacity)
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
        free(arr);
    }

    template<class T, size_t capacity>
    CStyleArray<T, capacity>& CStyleArray<T, capacity>::operator=(const CStyleArray<T, capacity>& array) {
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
    CStyleArray<T, capacity>&
    CStyleArray<T, capacity>::operator=(CStyleArray<T, capacity>&& array) noexcept {
        this->~CStyleArray();
        arr = array.arr;
        array.arr = nullptr;
        length = array.length;
        array.length = 0;
        return *this;
    }

    template<class T, size_t capacity>
    bool CStyleArray<T, capacity>::operator==(const CStyleArray<T, capacity>& array) const {
        if(length != array.length)
            return false;
        for(size_t i = 0; i < length; ++i)
            if(operator[](i) != array[i])
                return false;
        return true;
    }
    /*!
     * Return the sub array of current array. \from is included and \to is excluded.
     */
    template<class T, size_t capacity>
    CStyleArray<T, Dynamic> CStyleArray<T, capacity>::subArray(size_t from, size_t to) {
        Q_ASSERT(from < to && to <= length);
        const auto result_length = to - from;
        CStyleArray<T, Dynamic> result(result_length);
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < result_length; ++i, ++from)
                new (result.arr + i) T(arr[from]);
        else
            memcpy(result.arr, arr, result_length * sizeof(T));
    }
    /*!
     * Cut the array from index \from, \from is included.
     * The result is a array whose length and capacity are equal.
     */
    template<class T, size_t capacity>
    CStyleArray<T, Dynamic> CStyleArray<T, capacity>::cut(size_t from) {
        Q_ASSERT(from < length);
        auto result_length = length - from;
        CStyleArray<T, Dynamic> result(result_length);
        for(size_t i = 0; from < length; ++from, ++i)
            result.allocate(std::move(arr[from]), i);
        length = from;
        return result;
    }
    /*!
     * Low level api. Designed for performance.
     * Simply allocate a \T at position \index.
     * You must ensure position \index is usable or a memory leak will occur.
     */
    template<class T, size_t capacity>
    inline void CStyleArray<T, capacity>::allocate(const T& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(t);
        else
            *(arr + index) = t;
    }

    template<class T, size_t capacity>
    inline void CStyleArray<T, capacity>::allocate(T&& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(std::move(t));
        else
            *(arr + index) = t;
    }
    /*!
     * append() will add a element or array to the end of this array.
     * Wrap structure: append() <- grow() <- allocate()
     */
    template<class T, size_t capacity>
    void CStyleArray<T, capacity>::append(const T& t) {
        if(length == capacity)
            return;
        grow(t);
    }

    template<class T, size_t capacity>
    void CStyleArray<T, capacity>::append(T&& t) {
        if(length == capacity)
            return;
        grow(std::move(t));
    }

    template<class T, size_t capacity>
    void CStyleArray<T, capacity>::append(const CStyleArray<T, capacity>& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        t_length = new_length > capacity ? capacity - length : t_length;
        for(size_t i = 0; i < t_length; ++i, ++length)
            allocate(t[i], length);
    }

    template<class T, size_t capacity>
    void CStyleArray<T, capacity>::append(CStyleArray<T, capacity>&& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        const bool overflow = new_length > capacity;
        t_length = overflow ? capacity - length : t_length;

        memcpy(arr + length, t.arr, t_length * sizeof(T));
        length = overflow ? capacity : new_length;
        t.arr = nullptr;
        t.length = 0;
    }
    /*!
     * Low level api. Designed for performance.
     * Allocate a element at the end and increase the length.
     * This function can be used when you are sure the current capacity is enough.
     */
    template<class T, size_t capacity>
    inline void CStyleArray<T, capacity>::grow(const T& t) {
        Q_ASSERT(length < capacity);
        allocate(t, length++);
    }

    template<class T, size_t capacity>
    inline void CStyleArray<T, capacity>::grow(T&& t) {
        Q_ASSERT(length < capacity);
        allocate(std::move(t), length++);
    }
    /*!
     * Get the last element in the array and remove it from the array.
     */
    template<class T, size_t capacity>
    T CStyleArray<T, capacity>::cutLast() {
        Q_ASSERT(length > 0);
        --length;
        if(QTypeInfo<T>::isComplex)
            return T(std::move(arr[length]));
        else
            return arr[length];
    }

    template<class T, size_t capacity>
    void CStyleArray<T, capacity>::swap(CStyleArray<T, capacity>& array) noexcept {
        std::swap(arr, array.arr);
        std::swap(length, array.length);
    }
    ///////////////////////////////////////CStyleArray<T, Dynamic>//////////////////////////////////////////

    template<class T>
    inline CStyleArray<T, Dynamic>::CStyleArray()
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))), length(0) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not AbstractCStyleArray<T>::allocated. DO NOT try to visit them.
     */
    template<class T>
    inline CStyleArray<T, Dynamic>::CStyleArray(size_t length)
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))), length(length), capacity(length) {}

    template<class T>
    inline CStyleArray<T, Dynamic>::CStyleArray(const CStyleArray<T, Dynamic>& array)
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))), length(array.length) {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                new (arr + i) T(array.arr[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T>
    inline CStyleArray<T, Dynamic>::CStyleArray(CStyleArray<T, Dynamic>&& array) noexcept
            : arr(array.arr), length(array.length) {
        array.arr = nullptr;
        array.length = 0;
    }

    template<class T>
    CStyleArray<T, Dynamic>::~CStyleArray() {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
        free(arr);
    }

    template<class T>
    CStyleArray<T, Dynamic>& CStyleArray<T, Dynamic>::operator=(const CStyleArray<T, Dynamic>& array) {
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

    template<class T>
    CStyleArray<T, Dynamic>& CStyleArray<T, Dynamic>::operator=(CStyleArray<T, Dynamic>&& array) noexcept {
        this->~CStyleArray();
        arr = array.arr;
        array.arr = nullptr;
        length = array.length;
        array.length = 0;
        return *this;
    }

    template<class T>
    bool CStyleArray<T, Dynamic>::operator==(const CStyleArray<T, Dynamic>& array) const {
        if(length != array.length)
            return false;
        for(size_t i = 0; i < length; ++i)
            if(operator[](i) != array[i])
                return false;
        return true;
    }
    /*!
     * Return the sub array of current array. \from is included and \to is excluded.
     */
    template<class T>
    CStyleArray<T, Dynamic> CStyleArray<T, Dynamic>::subArray(size_t from, size_t to) {
        Q_ASSERT(from < to && to <= length);
        const auto result_length = to - from;
        CStyleArray<T, Dynamic> result(result_length);
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < result_length; ++i, ++from)
                new (result.arr + i) T(arr[from]);
        else
            memcpy(result.arr, arr, result_length * sizeof(T));
    }
    /*!
     * Cut the array from index \from, \from is included.
     * The result is a array whose length and capacity are equal.
     */
    template<class T>
    CStyleArray<T, Dynamic> CStyleArray<T, Dynamic>::cut(size_t from) {
        Q_ASSERT(from < length);
        auto result_length = length - from;
        CStyleArray<T, Dynamic> result(result_length);
        for(size_t i = 0; from < length; ++from, ++i)
            result.allocate(std::move(arr[from]), i);
        length = from;
        return result;
    }
    /*!
     * Low level api. Designed for performance.
     * Simply allocate a \T at position \index.
     * You must ensure position \index is usable or a memory leak will occur.
     */
    template<class T>
    inline void CStyleArray<T, Dynamic>::allocate(const T& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(t);
        else
            *(arr + index) = t;
    }

    template<class T>
    inline void CStyleArray<T, Dynamic>::allocate(T&& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(std::move(t));
        else
            *(arr + index) = t;
    }
    template<class T>
    void CStyleArray<T, Dynamic>::append(const T& t) {
        if(length == capacity)
            increase(capacity + 1);
        grow(t);
    }

    template<class T>
    void CStyleArray<T, Dynamic>::append(T&& t) {
        if(length == capacity)
            increase(capacity + 1);
        grow(std::move(t));
    }

    template<class T>
    void CStyleArray<T, Dynamic>::append(const CStyleArray& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        if(new_length > capacity)
            resize(new_length);
        for(size_t i = 0; i < t_length; ++i, ++length)
            allocate(t[i], length);
    }

    template<class T>
    void CStyleArray<T, Dynamic>::append(CStyleArray&& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        if(new_length > capacity)
            resize(new_length);
        memcpy(arr + length, t.arr, t_length * sizeof(T));
        length = new_length;
        t.arr = nullptr;
        t.length = 0;
    }
    /*!
     * Low level api. Designed for performance.
     * Allocate a element at the end and increase the length.
     * This function can be used when you are sure the current capacity is enough.
     */
    template<class T>
    inline void CStyleArray<T, Dynamic>::grow(const T& t) {
        Q_ASSERT(length < capacity);
        allocate(t, length++);
    }

    template<class T>
    inline void CStyleArray<T, Dynamic>::grow(T&& t) {
        Q_ASSERT(length < capacity);
        allocate(std::move(t), length++);
    }

    template<class T>
    void CStyleArray<T, Dynamic>::resize(size_t size) {
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

    template<class T>
    void CStyleArray<T, Dynamic>::squeeze() {
        arr = reinterpret_cast<T*>(realloc(arr, length * sizeof(T)));
        capacity = length;
    }
    /*!
     * Increase the capacity.
     * This function can be used when you are sure the new \size is larger than the old capacity.
     */
    template<class T>
    void CStyleArray<T, Dynamic>::increase(size_t size) {
        Q_ASSERT(size >= capacity);
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        capacity = size;
    }
    /*!
     * Decrease the capacity.
     * This function can be used when you are sure the new \size is shorter than the old capacity.
     */
    template<class T>
    void CStyleArray<T, Dynamic>::decrease(size_t size) {
        Q_ASSERT(size <= capacity);
        if(QTypeInfo<T>::isComplex) {
            for(size_t i = size; size < length; ++i)
                (arr + i)->~T();
        }
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        length = capacity = size;
    }
    /*!
     * Get the last element in the array and remove it from the array.
     */
    template<class T>
    T CStyleArray<T, Dynamic>::cutLast() {
        Q_ASSERT(length > 0);
        --length;
        if(QTypeInfo<T>::isComplex)
            return T(std::move(arr[length]));
        else
            return arr[length];
    }

    template<class T>
    void CStyleArray<T, Dynamic>::swap(CStyleArray& array) noexcept {
        std::swap(arr, array.arr);
        std::swap(length, array.length);
        std::swap(capacity, array.capacity);
    }
}

#endif