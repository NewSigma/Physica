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

namespace Physica::Utils {
    //////////////////////////////////////////CStyleArray<T, Length, Capacity>//////////////////////////////////////////
    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>::CStyleArray() : Base(Capacity) {
        for (size_t i = 0; i < Length; ++i)
            Base::allocate(T(), i);
    }

    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>::CStyleArray(std::initializer_list<T> list) : Base(Length) {
        static_assert(list.size() <= Capacity);
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            Base::allocate(*ite, i);
        for (; i < Length; ++i)
            Base::allocate(T(), i);
    }

    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>::CStyleArray(const CStyleArray<T, Length, Capacity>& array) : Base(array) {}

    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Length, Capacity>::CStyleArray(CStyleArray<T, Length, Capacity>&& array) noexcept
            : Base(std::move(array)) {}

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
    CStyleArray<T, Dynamic, Dynamic> CStyleArray<T, Length, Capacity>::subArray(size_t from, size_t to) {
        assert(from < to && to <= Length);
        const auto result_length = to - from;
        CStyleArray<T, Dynamic, Dynamic> result(result_length);
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < result_length; ++i, ++from)
                new (result.arr + i) T(arr[from]);
        else
            memcpy(result.arr, arr, result_length * sizeof(T));
        return result;
    }
    /*!
     * Cut the array from index \from, \from is included.
     * The result is a array whose length and capacity are equal.
     */
    template<class T, size_t Length, size_t Capacity>
    CStyleArray<T, Dynamic, Dynamic> CStyleArray<T, Length, Capacity>::cut(size_t from) {
        assert(from < Length);
        auto result_length = Length - from;
        CStyleArray<T, Dynamic, Dynamic> result(result_length);
        for(size_t i = 0; from < Length; ++from, ++i)
            result.allocate(std::move(Base::operator[](from)), i);
        return result;
    }
    ///////////////////////////////////////CStyleArray<T, Dynamic, Capacity>//////////////////////////////////////////
    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray() : Base(0, Capacity) {}
    /**
     * Low level api. Designed for performance.
     * Warning: The first \param length elements have not been allocated. DO NOT try to visit them.
     */
    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray(size_t length) : Base(length, Capacity) {}

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray(std::initializer_list<T> list) : Base(list.size(), Capacity) {
        static_assert(list.size() <= Capacity);
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            Base::allocate(*ite, i);
    }

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray(const CStyleArray<T, Dynamic, Capacity>& array)
            : Base(array) {}

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::CStyleArray(CStyleArray<T, Dynamic, Capacity>&& array) noexcept
            : Base(std::move(array)) {}

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>::~CStyleArray() {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
    }

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>& CStyleArray<T, Dynamic, Capacity>::operator=(const CStyleArray<T, Dynamic, Capacity>& array) {
        if(this != &array)
            Base::operator=(array);
        return *this;
    }

    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Capacity>& CStyleArray<T, Dynamic, Capacity>::operator=(CStyleArray<T, Dynamic, Capacity>&& array) noexcept {
        Base::operator=(std::move(array));
        return *this;
    }
    /**
     * Return the sub array of current array. \from is included and \to is excluded.
     */
    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Dynamic> CStyleArray<T, Dynamic, Capacity>::subArray(size_t from, size_t to) {
        assert(from < to && to <= length);
        const auto result_length = to - from;
        CStyleArray<T, Dynamic, Dynamic> result(result_length);
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < result_length; ++i, ++from)
                new (result.arr + i) T(arr[from]);
        else
            memcpy(result.arr, arr, result_length * sizeof(T));
        return result;
    }
    /**
     * Cut the array from index \from, \from is included.
     * The result is a array whose length and capacity are equal.
     */
    template<class T, size_t Capacity>
    CStyleArray<T, Dynamic, Dynamic> CStyleArray<T, Dynamic, Capacity>::cut(size_t from) {
        assert(from < length);
        auto result_length = length - from;
        CStyleArray<T, Dynamic, Dynamic> result(result_length);
        for(size_t i = 0; from < length; ++from, ++i)
            result.allocate(std::move(arr[from]), i);
        length = from;
        return result;
    }
    /**
     * append() will add a element or array to the end of this array.
     * Wrap structure: append() <- grow() <- allocate()
     */
    template<class T, size_t Capacity>
    inline void CStyleArray<T, Dynamic, Capacity>::append(const T& t) {
        assert(length < Capacity);
        Base::grow(t);
    }

    template<class T, size_t Capacity>
    inline void CStyleArray<T, Dynamic, Capacity>::append(T&& t) {
        assert(length < Capacity);
        Base::grow(std::move(t));
    }

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::append(const CStyleArray& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        t_length = new_length > Capacity ? Capacity - length : t_length;
        for(size_t i = 0; i < t_length; ++i, ++length)
            Base::allocate(t[i], length);
    }

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::append(CStyleArray&& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        const bool overflow = new_length > Capacity;
        t_length = overflow ? Capacity - length : t_length;

        memcpy(arr + length, t.arr, t_length * sizeof(T));
        length = overflow ? Capacity : new_length;
        t.arr = nullptr;
        t.length = 0;
    }

    template<class T, size_t Capacity>
    void CStyleArray<T, Dynamic, Capacity>::swap(CStyleArray<T, Dynamic, Capacity>& array) noexcept {
        Base::swap(array);
    }
    ///////////////////////////////////////CStyleArray<T, Dynamic, Dynamic>//////////////////////////////////////////
    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::CStyleArray() : Base(0, 0), capacity(0) {}

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::CStyleArray(size_t length_) : CStyleArray(length_, length_) {}
    /**
     * Low level api. Designed for performance.
     * Warning: The first \param length elements have not been allocated. DO NOT try to visit them.
     */
    template<class T>
    inline CStyleArray<T, Dynamic, Dynamic>::CStyleArray(size_t length_, size_t capacity_)
            : Base(length_, capacity_), capacity(capacity_) {
        assert(length_ <= capacity_);
    }

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::CStyleArray(std::initializer_list<T> list)
            : Base(list.size(), list.size()), capacity(list.size()) {
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            allocate(*ite, i);
    }

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::CStyleArray(const CStyleArray& array)
            : Base(array), capacity(array.capacity) {}

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::CStyleArray(CStyleArray<T, Dynamic, Dynamic>&& array) noexcept
            : Base(std::move(array)), capacity(array.capacity) {}

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>::~CStyleArray() {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
    }

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>& CStyleArray<T, Dynamic, Dynamic>::operator=(const CStyleArray<T, Dynamic, Dynamic>& array) {
        if(this != &array) {
            capacity = array.capacity;
            Base::operator=(array);
        }
        return *this;
    }

    template<class T>
    CStyleArray<T, Dynamic, Dynamic>& CStyleArray<T, Dynamic, Dynamic>::operator=(CStyleArray<T, Dynamic, Dynamic>&& array) noexcept {
        capacity = array.capacity;
        Base::operator=(std::move(array));
        return *this;
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(const T& t) {
        if(length == capacity)
            increase(capacity + 1);
        Base::grow(t);
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(T&& t) {
        if(length == capacity)
            increase(capacity + 1);
        Base::grow(std::move(t));
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(const CStyleArray& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        if(new_length > capacity)
            resize(new_length);
        for(size_t i = 0; i < t_length; ++i, ++length)
            Base::allocate(t[i], length);
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::append(CStyleArray&& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        if(new_length > capacity)
            resize(new_length);
        memcpy(arr + length, t.arr, t_length * sizeof(T));
        length = new_length;
        t.arr = nullptr;
        t.length = 0;
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::resize(size_t size) {
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
    void CStyleArray<T, Dynamic, Dynamic>::squeeze() {
        arr = reinterpret_cast<T*>(realloc(arr, length * sizeof(T)));
        capacity = length;
    }
    /*!
     * Increase the capacity.
     * This function can be used when you are sure the new \param size is larger than the old capacity.
     */
    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::increase(size_t size) {
        assert(size >= capacity);
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        capacity = size;
    }
    /*!
     * Decrease the capacity.
     * This function can be used when you are sure the new \param size is shorter than the old capacity.
     */
    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::decrease(size_t size) {
        assert(size <= capacity);
        if(QTypeInfo<T>::isComplex) {
            for(size_t i = size; i < length; ++i)
                (arr + i)->~T();
        }
        arr = reinterpret_cast<T*>(realloc(arr, size * sizeof(T)));
        length = capacity = size;
    }

    template<class T>
    void CStyleArray<T, Dynamic, Dynamic>::swap(CStyleArray& array) noexcept {
        Base::swap(array);
        std::swap(capacity, array.capacity);
    }
}
