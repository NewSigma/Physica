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
    //////////////////////////////////////////Array<T, Length, Capacity, Allocator>//////////////////////////////////////////
    template<class T, size_t Length, size_t Capacity, class Allocator>
    Array<T, Length, Capacity, Allocator>::Array() : Base(Capacity) {
        for (size_t i = 0; i < Length; ++i)
            alloc.construct(arr + i);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    Array<T, Length, Capacity, Allocator>::Array(size_t length_, ConstLValueReferenceType t) : Base(length_) {
        assert(length_ == Length);
        for (size_t i = 0; i < Length; ++i)
            alloc.construct(arr + i, t);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    Array<T, Length, Capacity, Allocator>::Array(std::initializer_list<T> list) : Base(Length) {
        assert(list.size() <= Capacity);
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            alloc.construct(arr + i, *ite);
        for (; i < Length; ++i)
            alloc.construct(arr + i);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    Array<T, Length, Capacity, Allocator>::Array(const Array<T, Length, Capacity, Allocator>& array) : Base(array) {}

    template<class T, size_t Length, size_t Capacity, class Allocator>
    Array<T, Length, Capacity, Allocator>::Array(Array<T, Length, Capacity, Allocator>&& array) noexcept
            : Base(std::move(array)) {}

    template<class T, size_t Length, size_t Capacity, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Length, Capacity, Allocator>::Array(const Array<T, OtherLength, OtherCapacity, Allocator>& array) : Base(Capacity) {
        static_assert(OtherLength == Dynamic, "OtherLength must equal to Length or be Dynamic.");
        assert(array.getLength() == Length);
        if constexpr (!std::is_trivial<T>::value)
            for (size_t i = 0; i < Length; ++i)
                alloc.allocate(arr + i, array[i]);
        else
            memcpy(arr, array.arr, Length * sizeof(T));
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Length, Capacity, Allocator>::Array(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept
            : Base(alloc.reallocate(array.arr, Length)) {
        static_assert(OtherLength == Dynamic, "OtherLength must equal to Length or be Dynamic.");
        assert(array.getLength() == Length);
        array.arr = nullptr;
        array.setLength(0);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Length, Capacity, Allocator>::Array(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array) : Base(Capacity) {
        assert(Length == array.getLength());
        for (size_t i = 0; i < Length; ++i)
            alloc.construct(arr + i, array[i]);
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Length, Capacity, Allocator>::Array(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept {
        assert(Length == array.getLength());
        for (size_t i = 0; i < Length; ++i)
            alloc.construct(arr + i, std::move(array[i]));
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    Array<T, Length, Capacity, Allocator>&
    Array<T, Length, Capacity, Allocator>::operator=(const Array<T, Length, Capacity, Allocator>& array) {
        if (this != &array) {
            if constexpr (!std::is_trivial<T>::value)
                for (size_t i = 0; i < Length; ++i)
                    this->operator[](i) = array[i];
            else
                memcpy(arr, array.arr, Length * sizeof(T));
        }
        return *this;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    Array<T, Length, Capacity, Allocator>&
    Array<T, Length, Capacity, Allocator>::operator=(Array<T, Length, Capacity, Allocator>&& array) noexcept {
        this->~Array();
        arr = array.arr;
        array.arr = nullptr;
        return *this;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Length, Capacity, Allocator>& Array<T, Length, Capacity, Allocator>::operator=(const Array<T, OtherLength, OtherCapacity, Allocator>& array) {
        static_assert(OtherLength == Dynamic, "OtherLength must equal to Length or be Dynamic.");
        assert(array.getLength() == Length);
        if constexpr (!std::is_trivial<T>::value)
            for (size_t i = 0; i < Length; ++i)
                this->operator[](i) = array[i];
        else
            memcpy(arr, array.arr, Length * sizeof(T));
        return *this;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Length, Capacity, Allocator>& Array<T, Length, Capacity, Allocator>::operator=(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept {
        static_assert(OtherLength == Dynamic, "OtherLength must equal to Length or be Dynamic.");
        assert(array.getLength() == Length);
        this->~Array();
        arr = alloc.reallocate(array.arr, Length);
        array.arr = nullptr;
        array.setLength(0);
        return *this;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Length, Capacity, Allocator>& Array<T, Length, Capacity, Allocator>::operator=(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array) {
        assert(Length == array.getLength());
        for (size_t i = 0; i < Length; ++i)
            this->operator[](i) = array[i];
        return *this;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Length, Capacity, Allocator>& Array<T, Length, Capacity, Allocator>::operator=(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept {
        assert(Length == array.getLength());
        for (size_t i = 0; i < Length; ++i)
            this->operator[](i) = std::move(array[i]);
        return *this;
    }
    /*!
     * Return the sub array of current array. \from is included and \to is excluded.
     */
    template<class T, size_t Length, size_t Capacity, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator> Array<T, Length, Capacity, Allocator>::subArray(size_t from, size_t to) {
        assert(from < to && to <= Length);
        const auto result_length = to - from;
        Array<T, Dynamic, Dynamic, Allocator> result(result_length);
        if constexpr(!std::is_trivial<T>::value)
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
    template<class T, size_t Length, size_t Capacity, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator> Array<T, Length, Capacity, Allocator>::cut(size_t from) {
        assert(from < Length);
        auto result_length = Length - from;
        Array<T, Dynamic, Dynamic, Allocator> result(result_length);
        for(size_t i = 0; from < Length; ++from, ++i)
            result.allocate(std::move(Base::operator[](from)), i);
        return result;
    }

    template<class T, size_t Length, size_t Capacity, class Allocator>
    void Array<T, Length, Capacity, Allocator>::resize(size_t size, const T& t) {
        assert(size == getLength());
        for (size_t i = 0; i < size; ++i)
            (*this)[i] = t;
    }
    ///////////////////////////////////////Array<T, Dynamic, Capacity>//////////////////////////////////////////
    template<class T, size_t Capacity, class Allocator>
    Array<T, Dynamic, Capacity, Allocator>::Array() : Base(Capacity) {}

    template<class T, size_t Capacity, class Allocator>
    Array<T, Dynamic, Capacity, Allocator>::Array(size_t length_, ConstLValueReferenceType t) : Base(length_, Capacity) {
        assert(length_ < Capacity);
        for (size_t i = 0; i < length_; ++i)
            alloc.construct(arr + i, t);
    }

    template<class T, size_t Capacity, class Allocator>
    Array<T, Dynamic, Capacity, Allocator>::Array(std::initializer_list<T> list) : Base(Capacity) {
        constexpr auto length = list.size();
        static_assert(length <= Capacity);
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            alloc.construct(arr + i, *ite);
        Base::setLength(length);
    }

    template<class T, size_t Capacity, class Allocator>
    Array<T, Dynamic, Capacity, Allocator>::Array(const Array<T, Dynamic, Capacity, Allocator>& array)
            : Base(array) {}

    template<class T, size_t Capacity, class Allocator>
    Array<T, Dynamic, Capacity, Allocator>::Array(Array<T, Dynamic, Capacity, Allocator>&& array) noexcept
            : Base(std::move(array)) {}

    template<class T, size_t Capacity, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Capacity, Allocator>::Array(const Array<T, OtherLength, OtherCapacity, Allocator>& array) : Base(array.getLength(), Capacity) {
        assert(array.getLength() <= Capacity);
        if constexpr (!std::is_trivial<T>::value)
            for (size_t i = 0; i < length; ++i)
                alloc.construct(arr + i, array[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T, size_t Capacity, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Capacity, Allocator>::Array(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept
            : Base(array.getLength(), alloc.reallocate(array.arr, Capacity)) {
        assert(array.getLength() <= Capacity);
        array.arr = nullptr;
        array.setLength(0);
    }

    template<class T, size_t Capacity, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Capacity, Allocator>::Array(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array) : Base(array.getLength(), Capacity) {
        for (size_t i = 0; i < length; ++i)
            alloc.construct(arr + i, array[i]);
    }

    template<class T, size_t Capacity, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Capacity, Allocator>::Array(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept : Base(array.getLength(), Capacity) {
        for (size_t i = 0; i < length; ++i)
            alloc.construct(arr + i, std::move(array[i]));
    }

    template<class T, size_t Capacity, class Allocator>
    Array<T, Dynamic, Capacity, Allocator>&
    Array<T, Dynamic, Capacity, Allocator>::operator=(const Array<T, Dynamic, Capacity, Allocator>& array) {
        if (this != &array) {
            Base::clear();
            Base::setLength(array.getLength());
            if constexpr (!std::is_trivial<T>::value)
                for (size_t i = 0; i < length; ++i)
                    this->operator[](i) = array[i];
            else
                memcpy(arr, array.arr, length * sizeof(T));
        }
        return *this;
    }

    template<class T, size_t Capacity, class Allocator>
    Array<T, Dynamic, Capacity, Allocator>&
    Array<T, Dynamic, Capacity, Allocator>::operator=(Array<T, Dynamic, Capacity, Allocator>&& array) noexcept {
        this->~Array();
        arr = array.arr;
        array.arr = nullptr;
        length = array.length;
        array.length = 0;
        return *this;
    }

    template<class T, size_t Capacity, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Capacity, Allocator>&
    Array<T, Dynamic, Capacity, Allocator>::operator=(const Array<T, OtherLength, OtherCapacity, Allocator>& array) {
        assert(array.getLength() <= Capacity);
        Base::clear();
        Base::setLength(array.getLength());
        if constexpr (!std::is_trivial<T>::value)
            for (size_t i = 0; i < length; ++i)
                alloc.construct(arr + i, array[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
        return *this;
    }

    template<class T, size_t Capacity, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Capacity, Allocator>&
    Array<T, Dynamic, Capacity, Allocator>::operator=(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept {
        assert(array.getLength() <= Capacity);
        this->~Array();
        arr = alloc.reallocate(array.arr, Capacity);
        Base::setLength(array.getLength());
        array.arr = nullptr;
        array.setLength(0);
        return *this;
    }

    template<class T, size_t Capacity, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Capacity, Allocator>&
    Array<T, Dynamic, Capacity, Allocator>::operator=(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array) {
        assert(array.getLength() <= Capacity);
        Base::clear();
        for (size_t i = 0; i < array.getLength(); ++i)
            alloc.construct(arr + i, array[i]);
        setLength(array.getLength());
        return *this;
    }

    template<class T, size_t Capacity, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Capacity, Allocator>&
    Array<T, Dynamic, Capacity, Allocator>::operator=(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept {
        assert(array.getLength() <= Capacity);
        Base::clear();
        for (size_t i = 0; i < array.getLength(); ++i)
            alloc.construct(arr + i, std::move(array[i]));
        setLength(array.getLength());
        return *this;
    }
    /**
     * Return the sub array of current array. \from is included and \to is excluded.
     */
    template<class T, size_t Capacity, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator> Array<T, Dynamic, Capacity, Allocator>::subArray(size_t from, size_t to) {
        assert(from < to && to <= length);
        const auto result_length = to - from;
        Array<T, Dynamic, Dynamic, Allocator> result(result_length);
        if(!std::is_trivial<T>::value)
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
    template<class T, size_t Capacity, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator> Array<T, Dynamic, Capacity, Allocator>::cut(size_t from) {
        assert(from < length);
        auto result_length = length - from;
        Array<T, Dynamic, Dynamic, Allocator> result(result_length);
        for(size_t i = 0; from < length; ++from, ++i)
            result.allocate(std::move(arr[from]), i);
        length = from;
        return result;
    }
    /**
     * append() will add a element or array to the end of this array.
     * Wrap structure: append() <- grow() <- allocate()
     */
    template<class T, size_t Capacity, class Allocator>
    inline void Array<T, Dynamic, Capacity, Allocator>::append(ConstLValueReferenceType t) {
        assert(length < Capacity);
        Base::grow(t);
    }

    template<class T, size_t Capacity, class Allocator>
    inline void Array<T, Dynamic, Capacity, Allocator>::append(RValueReferenceType t) {
        assert(length < Capacity);
        Base::grow(std::move(t));
    }

    template<class T, size_t Capacity, class Allocator>
    void Array<T, Dynamic, Capacity, Allocator>::append(const Array& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        t_length = new_length > Capacity ? Capacity - length : t_length;
        for(size_t i = 0; i < t_length; ++i, ++length)
            alloc.construct(arr + length, t[i]);
    }

    template<class T, size_t Capacity, class Allocator>
    void Array<T, Dynamic, Capacity, Allocator>::append(Array&& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        const bool overflow = new_length > Capacity;
        t_length = overflow ? Capacity - length : t_length;

        memcpy(arr + length, t.arr, t_length * sizeof(T));
        length = overflow ? Capacity : new_length;
        t.arr = nullptr;
        t.length = 0;
    }
    /**
     * For the convenience of implementing templates.
     */
    template<class T, size_t Capacity, class Allocator>
    void Array<T, Dynamic, Capacity, Allocator>::reserve(size_t size) {
        assert(size == Capacity);
    }

    template<class T, size_t Capacity, class Allocator>
    void Array<T, Dynamic, Capacity, Allocator>::swap(Array<T, Dynamic, Capacity, Allocator>& array) noexcept {
        Base::swap(array);
    }
    ///////////////////////////////////////Array<T, Dynamic, Dynamic, Allocator>//////////////////////////////////////////
    template<class T, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator>::Array() : Base(0), capacity(0) {}

    template<class T, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator>::Array(size_t length_, ConstLValueReferenceType t) : Base(length_, length_), capacity(length_) {
        for (size_t i = 0; i < length_; ++i)
            alloc.construct(arr + i, t);
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator>::Array(std::initializer_list<T> list)
            : Base(list.size()), capacity(list.size()) {
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            alloc.construct(arr + i, *ite);
        Base::setLength(list.size());
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator>::Array(const Array& array)
            : Base(array), capacity(array.capacity) {}

    template<class T, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator>::Array(Array<T, Dynamic, Dynamic, Allocator>&& array) noexcept
            : Base(std::move(array)), capacity(array.capacity) {}

    template<class T, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Dynamic, Allocator>::Array(const Array<T, OtherLength, OtherCapacity, Allocator>& array)
            : Base(array.getLength(), array.getCapacity()), capacity(array.getCapacity()) {
        if constexpr (!std::is_trivial<T>::value)
            for (size_t i = 0; i < length; ++i)
                alloc.construct(arr + i, array[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class T, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Dynamic, Allocator>::Array(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept
            : Base(array.getLength(), array.arr), capacity(array.getCapacity()) {
        array.arr = nullptr;
        array.setLength(0);
    }

    template<class T, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Dynamic, Allocator>::Array(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array)
            : Base(array.getLength(), array.getLength()), capacity(array.getLength()) {
        for (size_t i = 0; i < array.getLength(); ++i)
            alloc.construct(arr + i, array[i]);
    }

    template<class T, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Dynamic, Allocator>::Array(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept
            : Base(array.getLength(), array.getLength()), capacity(array.getLength()) {
        for (size_t i = 0; i < array.getLength(); ++i)
            alloc.construct(arr + i, std::move(array[i]));
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator>&
    Array<T, Dynamic, Dynamic, Allocator>::operator=(const Array<T, Dynamic, Dynamic, Allocator>& array) {
        if (this != &array) {
            this->~Array();
            capacity = array.capacity;
            arr = reinterpret_cast<T*>(malloc(capacity * sizeof(T)));
            length = array.length;
            if constexpr (!std::is_trivial<T>::value)
                for (size_t i = 0; i < length; ++i)
                    this->operator[](i) = array[i];
            else
                memcpy(arr, array.arr, length * sizeof(T));
        }
        return *this;
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Dynamic, Allocator>&
    Array<T, Dynamic, Dynamic, Allocator>::operator=(Array<T, Dynamic, Dynamic, Allocator>&& array) noexcept {
        this->~Array();
        arr = array.arr;
        array.arr = nullptr;
        length = array.length;
        array.length = 0;
        capacity = array.capacity;
        return *this;
    }

    template<class T, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Dynamic, Allocator>&
    Array<T, Dynamic, Dynamic, Allocator>::operator=(const Array<T, OtherLength, OtherCapacity, Allocator>& array) {
        Base::clear();
        Base::setLength(array.getLength());
        capacity = array.getCapacity();
        arr = alloc.reallocate(arr, capacity);
        if constexpr (!std::is_trivial<T>::value)
            for (size_t i = 0; i < length; ++i)
                alloc.construct(arr + i, array[i]);
        else
            memcpy(arr, array.arr, length * sizeof(T));
        return *this;
    }

    template<class T, class Allocator>
    template<size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Dynamic, Allocator>&
    Array<T, Dynamic, Dynamic, Allocator>::operator=(Array<T, OtherLength, OtherCapacity, Allocator>&& array) noexcept {
        this->~Array();
        length = array.getLength();
        capacity = array.getCapacity();
        arr = array.arr;
        return *this;
    }

    template<class T, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Dynamic, Allocator>&
    Array<T, Dynamic, Dynamic, Allocator>::operator=(const Array<OtherT, OtherLength, OtherCapacity, Allocator>& array) {
        Base::clear();
        setLength(array.getLength());
        arr = alloc.reallocate(arr, length);
        for (size_t i = 0; i < length; ++i)
            alloc.construct(arr + i, array[i]);
        return *this;
    }

    template<class T, class Allocator>
    template<class OtherT, size_t OtherLength, size_t OtherCapacity>
    Array<T, Dynamic, Dynamic, Allocator>&
    Array<T, Dynamic, Dynamic, Allocator>::operator=(Array<OtherT, OtherLength, OtherCapacity, Allocator>&& array) noexcept {
        Base::clear();
        setLength(array.getLength());
        arr = alloc.reallocate(arr, length);
        for (size_t i = 0; i < length; ++i)
            alloc.construct(arr + i, std::move(array[i]));
        return *this;
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::append(ConstLValueReferenceType t) {
        if(length == capacity)
            increase(capacity + 1);
        Base::grow(t);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::append(RValueReferenceType t) {
        if(length == capacity)
            increase(capacity + 1);
        Base::grow(std::move(t));
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::append(const Array& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        if(new_length > capacity)
            reserve(new_length);
        for(size_t i = 0; i < t_length; ++i, ++length)
            alloc.construct(arr + length, t[i]);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::append(Array&& t) {
        const auto t_length = t.length;
        const auto new_length = length + t_length;
        if(new_length > capacity)
            reserve(new_length);
        memcpy(arr + length, t.arr, t_length * sizeof(T));
        length = new_length;
        t.arr = nullptr;
        t.length = 0;
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::reserve(size_t size) {
        assert (size > getCapacity());
        arr = alloc.reallocate(arr, size);
        capacity = size;
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::resize(size_t size) {
        if constexpr (!std::is_trivial<T>::value) {
            if (length > size) {
                for (size_t i = size; i < length; ++i)
                    (arr + i)->~T();
                length = size;
            }
        }
        if (capacity < size)
            reserve(size);
        for (; length < size; ++length)
            alloc.construct(arr + length);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::resize(size_t size, const T& t) {
        if constexpr (!std::is_trivial<T>::value) {
            if (length > size) {
                for (size_t i = size; i < length; ++i)
                    (arr + i)->~T();
                length = size;
            }
        }
        if (capacity < size)
            reserve(size);
        for (; length < size; ++length)
            alloc.construct(arr + length, t);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::squeeze() {
        arr = alloc.reallocate(arr, length);
        capacity = length;
    }
    /*!
     * Increase the capacity.
     * This function can be used when you are sure the new \param size is larger than the old capacity.
     */
    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::increase(size_t size) {
        assert(size >= capacity);
        arr = alloc.reallocate(arr, size);
        capacity = size;
    }
    /*!
     * Decrease the capacity.
     * This function can be used when you are sure the new \param size is shorter than the old capacity.
     */
    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::decrease(size_t size) {
        assert(size <= capacity);
        if(!std::is_trivial<T>::value) {
            for(size_t i = size; i < length; ++i)
                (arr + i)->~T();
        }
        arr = alloc.reallocate(arr, size);
        length = capacity = size;
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Dynamic, Allocator>::swap(Array& array) noexcept {
        Base::swap(array);
        std::swap(capacity, array.capacity);
    }
}
