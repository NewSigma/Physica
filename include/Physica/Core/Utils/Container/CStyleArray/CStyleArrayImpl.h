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
#ifndef PHYSICA_CSTYLEARRAYIMPL_H
#define PHYSICA_CSTYLEARRAYIMPL_H

#include <cstring>

namespace Physica::Core {
    //////////////////////////////////////////CStyleArray<T, capacity>//////////////////////////////////////////
    template<class T, size_t capacity>
    inline CStyleArray<T, capacity>::CStyleArray()
            : AbstractCStyleArray<T>(reinterpret_cast<T*>(malloc(capacity * sizeof(T))), 0) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not AbstractCStyleArray<T>::allocated. DO NOT try to visit them.
     */
    template<class T, size_t capacity>
    inline CStyleArray<T, capacity>::CStyleArray(size_t length)
            : AbstractCStyleArray<T>(reinterpret_cast<T*>(malloc(capacity * sizeof(T))), length) {}

    template<class T, size_t capacity>
    CStyleArray<T, capacity>::CStyleArray(std::initializer_list<T> list) : AbstractCStyleArray<T>(list) {
        Q_UNUSED(capacity)
        static_assert(list.size() <= capacity);
    }

    template<class T, size_t capacity>
    inline CStyleArray<T, capacity>::CStyleArray(const CStyleArray<T, capacity>& array)
            : AbstractCStyleArray<T>(reinterpret_cast<T*>(malloc(capacity * sizeof(T))), array.length) {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                new (arr + i) T(array.arr[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T, size_t capacity>
    inline CStyleArray<T, capacity>::CStyleArray(CStyleArray<T, capacity>&& array) noexcept
            : AbstractCStyleArray<T>(static_cast<Base&&>(array)) {}

    template<class T, size_t capacity>
    CStyleArray<T, capacity>& CStyleArray<T, capacity>::operator=(const CStyleArray<T, capacity>& array) {
        if(this != &array) {
            this->~CStyleArray();

            arr = reinterpret_cast<T*>(malloc(capacity * sizeof(T)));
            if(QTypeInfo<T>::isComplex)
                for(size_t i = 0; i < length; ++i)
                    new (arr + i) T(array[i]);
            else
                memcpy(arr, array.arr, length * sizeof(T));
        }
        return *this;
    }

    template<class T, size_t capacity>
    CStyleArray<T, capacity>&
    CStyleArray<T, capacity>::operator=(CStyleArray<T, capacity>&& array) noexcept {
        Base::operator=(static_cast<Base&&>(array));
        return *this;
    }
    /*!
     * Return the sub array of current array. \from is included and \to is excluded.
     */
    template<class T, size_t capacity>
    CStyleArray<T, Dynamic> CStyleArray<T, capacity>::subArray(size_t from, size_t to) {
        Q_UNUSED(capacity)
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
        Q_UNUSED(capacity)
        Q_ASSERT(from < length);
        auto result_length = length - from;
        CStyleArray<T, Dynamic> result(result_length);
        for(size_t i = 0; from < length; ++from, ++i)
            result.allocate(std::move(Base::operator[](from)), i);
        return result;
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
            Base::allocate(t[i], length);
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
        Base::allocate(t, length++);
    }

    template<class T, size_t capacity>
    inline void CStyleArray<T, capacity>::grow(T&& t) {
        Q_ASSERT(length < capacity);
        Base::allocate(std::move(t), length++);
    }

    template<class T, size_t capacity>
    void CStyleArray<T, capacity>::removeAt(size_t index) {
        Q_UNUSED(capacity)
        Q_ASSERT(index < length);
        if(QTypeInfo<T>::isComplex)
            (arr + index)->~T();
        --length;
        memmove(arr + index, arr + index + 1, (length - index) * sizeof(T));
    }
    ///////////////////////////////////////CStyleArray<T, Dynamic>//////////////////////////////////////////
    template<class T>
    inline CStyleArray<T, Dynamic>::CStyleArray()
            : AbstractCStyleArray<T>(nullptr, 0), capacity(0) {}
    /*!
     * Low level api. Designed for performance.
     * Warning: The first \length elements have not AbstractCStyleArray<T>::allocated. DO NOT try to visit them.
     */
    template<class T>
    inline CStyleArray<T, Dynamic>::CStyleArray(size_t length)
            : AbstractCStyleArray<T>(reinterpret_cast<T*>(malloc(length * sizeof(T))), length), capacity(length) {}

    template<class T>
    CStyleArray<T, Dynamic>::CStyleArray(std::initializer_list<T> list)
            : AbstractCStyleArray<T>(list), capacity(list.size()) {}

    template<class T>
    inline CStyleArray<T, Dynamic>::CStyleArray(const CStyleArray<T, Dynamic>& array)
            : AbstractCStyleArray<T>(reinterpret_cast<T*>(malloc(array.capacity * sizeof(T))), array.length)
            , capacity(array.capacity) {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                new (arr + i) T(array.arr[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T>
    inline CStyleArray<T, Dynamic>::CStyleArray(CStyleArray<T, Dynamic>&& array) noexcept
            : AbstractCStyleArray<T>(static_cast<Base&&>(array)), capacity(array.capacity) {}

    template<class T>
    CStyleArray<T, Dynamic>& CStyleArray<T, Dynamic>::operator=(const CStyleArray<T, Dynamic>& array) {
        if(this != &array) {
            this->~CStyleArray();

            length = array.length;
            capacity = array.capacity;
            arr = reinterpret_cast<T*>(malloc(capacity * sizeof(T)));
            if(QTypeInfo<T>::isComplex)
                for(size_t i = 0; i < length; ++i)
                    new (arr + i) T(array[i]);
            else
                memcpy(arr, array.arr, length * sizeof(T));
        }
        return *this;
    }

    template<class T>
    CStyleArray<T, Dynamic>& CStyleArray<T, Dynamic>::operator=(CStyleArray<T, Dynamic>&& array) noexcept {
        Base::operator=(static_cast<Base&&>(array));
        capacity = array.capacity;
        return *this;
    }

    template<class T>
    bool CStyleArray<T, Dynamic>::operator==(const CStyleArray<T, Dynamic>& array) const {
        return Base::operator==(array) && capacity == array.capacity;
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
            Base::allocate(t[i], length);
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
        Base::allocate(t, length++);
    }

    template<class T>
    inline void CStyleArray<T, Dynamic>::grow(T&& t) {
        Q_ASSERT(length < capacity);
        Base::allocate(std::move(t), length++);
    }

    template<class T>
    void CStyleArray<T, Dynamic>::resize(size_t size) {
        if(QTypeInfo<T>::isComplex) {
            if(length > size) {
                for(size_t i = size; i < length; ++i)
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
            for(size_t i = size; i < length; ++i)
                (arr + i)->~T();
        }
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        length = capacity = size;
    }

    template<class T>
    void CStyleArray<T, Dynamic>::removeAt(size_t index) {
        Q_ASSERT(index < length);
        if(QTypeInfo<T>::isComplex)
            (arr + index)->~T();
        --length;
        memmove(arr + index, arr + index + 1, (length - index) * sizeof(T));
    }

    template<class T>
    void CStyleArray<T, Dynamic>::swap(CStyleArray<T, Dynamic>& array) noexcept {
        Base::swap(array);
        std::swap(capacity, array.capacity);
    }
}

#endif