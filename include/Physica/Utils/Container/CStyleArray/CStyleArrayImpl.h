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
#include <cstring>
#include <qtypeinfo.h>

namespace Physica::Utils {
    //////////////////////////////////////////CStyleArray<T, Length, Capacity>//////////////////////////////////////////
    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>::CStyleArray() : Base(Capacity) {
        for (size_t i = 0; i < Length; ++i)
            allocate(T(), i);
    }

    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>::CStyleArray(std::initializer_list<T> list) : Base(Length) {
        static_assert(list.size() <= capacity);
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            allocate(*ite, i);
        for (; i < Length; ++i)
            allocate(T(), i);
    }

    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>::CStyleArray(const CStyleArray<T, Length, Capacity>& array) : Base(array) {}

    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>::CStyleArray(CStyleArray<T, Length, Capacity>&& array) noexcept
            : Base>(static_cast<Base&&>(array)) {}

    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>::~CStyleArray() {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < Length; ++i)
                (arr + i)->~T();
    }

    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>&
    CStyleArray<T, Length, Capacity>::operator=(const CStyleArray<T, Length, Capacity>& array) {
        Base::operator=(array);
        return *this;
    }

    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>&
    CStyleArray<T, Length, Capacity>::operator=(CStyleArray<T, Length, Capacity>&& array) noexcept {
        Base::operator=(std::move(array));
        return *this;
    }
    /*!
     * Return the sub array of current array. \from is included and \to is excluded.
     */
    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity> CStyleArray<T, Length, Capacity>::subArray(size_t from, size_t to) {
        Q_UNUSED(capacity)
        assert(from < to && to <= length);
        const auto result_length = to - from;
        CStyleArray<T, Dynamic, Capacity> result(result_length);
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
    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity> CStyleArray<T, Length, Capacity>::cut(size_t from) {
        Q_UNUSED(capacity)
        assert(from < length);
        auto result_length = length - from;
        CStyleArray<T, Dynamic, Capacity> result(result_length);
        for(size_t i = 0; from < length; ++from, ++i)
            result.allocate(std::move(Base::operator[](from)), i);
        return result;
    }
    ///////////////////////////////////////CStyleArray<T, Dynamic, Capacity>//////////////////////////////////////////
    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray() : Base(Capacity), length(0) {}
    /**
     * Low level api. Designed for performance.
     * Warning: The first \param length elements have not been allocated. DO NOT try to visit them.
     */
    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray(size_t length) : Base(Capacity), length(length) {}

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray(std::initializer_list<T> list) : length(list.size()) {
        static_assert(list.size() <= capacity);
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            allocate(*ite, i);
    }

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray(const CStyleArray<T, Dynamic, Capacity>& array)
            : Base(array), length(array.length) {}

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray(CStyleArray<T, Dynamic, Capacity>&& array) noexcept
            : Base(static_cast<Base&&>(array)), length(array.length) {}

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::~CStyleArray() {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
    }

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>& CStyleArray<T, Dynamic, Capacity>::operator=(const CStyleArray<T, Dynamic, Capacity>& array) {
        if(this != &array) {
            length = array.length;
            Base::operator=(array);
        }
        return *this;
    }

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>& CStyleArray<T, Dynamic, Capacity>::operator=(CStyleArray<T, Dynamic, Capacity>&& array) noexcept {
        length = array.length;
        Base::operator=(std::move(array));
        return *this;
    }

    template<class T, size_t Capacity>
    bool CStyleArray<T, Dynamic, Capacity>::operator==(const CStyleArray<T, Dynamic, Capacity>& array) const {
        return Base::operator==(array) && capacity == array.capacity;
    }
    /**
     * Get the last element in the array and remove it from the array.
     */
    template<class T, size_t Capacity>
    T CStyleArray<T, Dynamic, Capacity>::cutLast() {
        assert(length > 0);
        --length;
        if(QTypeInfo<T>::isComplex)
            return T(std::move(arr[length]));
        else
            return arr[length];
    }
    /**
     * Return the sub array of current array. \from is included and \to is excluded.
     */
    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity> CStyleArray<T, Dynamic, Capacity>::subArray(size_t from, size_t to) {
        assert(from < to && to <= length);
        const auto result_length = to - from;
        CStyleArray<T, Dynamic, Capacity> result(result_length);
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
    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity> CStyleArray<T, Dynamic, Capacity>::cut(size_t from) {
        assert(from < length);
        auto result_length = length - from;
        CStyleArray<T, Dynamic, Capacity> result(result_length);
        for(size_t i = 0; from < length; ++from, ++i)
            result.allocate(std::move(arr[from]), i);
        length = from;
        return result;
    }
    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::append(const T& t) {
        if(length == capacity)
            increase(capacity + 1);
        grow(t);
    }

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::append(T&& t) {
        if(length == capacity)
            increase(capacity + 1);
        grow(std::move(t));
    }

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::append(const CStyleArray& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        if(new_length > capacity)
            resize(new_length);
        for(size_t i = 0; i < t_length; ++i, ++length)
            Base::allocate(t[i], length);
    }

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::append(CStyleArray&& t) {
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
    template<class T, size_t Capacity>
    inline void CStyleArray<T, Dynamic, Capacity>::grow(const T& t) {
        assert(length < capacity);
        Base::allocate(t, length++);
    }

    template<class T, size_t Capacity>
    inline void CStyleArray<T, Dynamic, Capacity>::grow(T&& t) {
        assert(length < capacity);
        Base::allocate(std::move(t), length++);
    }

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::resize(size_t size) {
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

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::squeeze() {
        arr = reinterpret_cast<T*>(realloc(arr, length * sizeof(T)));
        capacity = length;
    }
    /*!
     * Increase the capacity.
     * This function can be used when you are sure the new \size is larger than the old capacity.
     */
    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::increase(size_t size) {
        assert(size >= capacity);
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        capacity = size;
    }
    /*!
     * Decrease the capacity.
     * This function can be used when you are sure the new \size is shorter than the old capacity.
     */
    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::decrease(size_t size) {
        assert(size <= capacity);
        if(QTypeInfo<T>::isComplex) {
            for(size_t i = size; i < length; ++i)
                (arr + i)->~T();
        }
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        length = capacity = size;
    }

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::removeAt(size_t index) {
        assert(index < length);
        if(QTypeInfo<T>::isComplex)
            (arr + index)->~T();
        --length;
        memmove(arr + index, arr + index + 1, (length - index) * sizeof(T));
    }

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::swap(CStyleArray<T, Dynamic, Capacity>& array) noexcept {
        Base::swap(array);
        std::swap(capacity, array.capacity);
    }
    ///////////////////////////////////////CStyleArray<T, Dynamic, Dynamic>//////////////////////////////////////////
    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::CStyleArray() : Base(0), length(0), capacity(0) {}
    /**
     * Low level api. Designed for performance.
     * Warning: The first \param length elements have not been allocated. DO NOT try to visit them.
     */
    template<class T>
    inline CStyleArray<T, Dynamic, Dynamic>::CStyleArray(size_t length)
            : Base(length), length(length), capacity(length) {}

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::CStyleArray(std::initializer_list<T> list) : length(list.size()), capacity(list.size()) {
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            allocate(*ite, i);
    }

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::CStyleArray(const CStyleArray& array)
            : Base(array), length(array.length), capacity(array.capacity) {}

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::CStyleArray(CStyleArray<T, Dynamic, Dynamic>&& array)
            : Base(std::move(array)), length(array.length), capacity(array.capacity) {}

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::~CStyleArray() {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
    }

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>& CStyleArray<T, Dynamic, Dynamic>::operator=(const CStyleArray<T, Dynamic, Dynamic>& array) {
        if(this != &array) {
            length = array.length;
            capacity = array.capacity;
            Base::operator=(array);
        }
        return *this;
    }

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>& CStyleArray<T, Dynamic, Dynamic>::operator=(CStyleArray<T, Dynamic, Dynamic>&& array) noexcept {
        length = array.length;
        capacity = array.capacity;
        Base::operator=(std::move(array));
        return *this;
    }
    /**
     * Get the last element in the array and remove it from the array.
     */
    template<class T>
    T CStyleArray<T, Dynamic, Dynamic>::cutLast() {
        assert(length > 0);
        --length;
        if(QTypeInfo<T>::isComplex)
            return T(std::move(arr[length]));
        else
            return arr[length];
    }
    /**
     * append() will add a element or array to the end of this array.
     * Wrap structure: append() <- grow() <- allocate()
     */
    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(const T& t) {
        if(length == capacity)
            return;
        grow(t);
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(T&& t) {
        if(length == capacity)
            return;
        grow(std::move(t));
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(const CStyleArray<T, Dynamic, Dynamic>& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        t_length = new_length > capacity ? capacity - length : t_length;
        for(size_t i = 0; i < t_length; ++i, ++length)
            Base::allocate(t[i], length);
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(CStyleArray<T, Dynamic, Dynamic>&& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        const bool overflow = new_length > capacity;
        t_length = overflow ? capacity - length : t_length;

        memcpy(arr + length, t.arr, t_length * sizeof(T));
        length = overflow ? capacity : new_length;
        t.arr = nullptr;
        t.length = 0;
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::swap(CStyleArray<T, Dynamic, Dynamic>& array) noexcept {
        Base::swap(array);
        std::swap(length, array.length);
        std::swap(capacity, array.capacity);
    }
}
