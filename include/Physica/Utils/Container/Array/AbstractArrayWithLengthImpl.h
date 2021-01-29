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

#include <qtypeinfo.h>

namespace Physica::Utils::Intenal {
    template<class Derived>
    AbstractArrayWithLength<Derived>::AbstractArrayWithLength(size_t capacity)
        : Base(capacity), length(0) {}

    template<class Derived>
    AbstractArrayWithLength<Derived>::AbstractArrayWithLength(size_t length_, size_t capacity)
        : Base(capacity), length(length_) {}

    template<class Derived>
    AbstractArrayWithLength<Derived>::AbstractArrayWithLength(
        const AbstractArrayWithLength& array) : Base(array), length(array.length) {}
    
    template<class Derived>
    AbstractArrayWithLength<Derived>::AbstractArrayWithLength(
        AbstractArrayWithLength&& array) noexcept : Base(std::move(array)), length(array.length) {}

    template<class Derived>
    AbstractArrayWithLength<Derived>::~AbstractArrayWithLength() {
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                (arr + i)->~T();
    }

    template<class Derived>
    AbstractArrayWithLength<Derived>&
    AbstractArrayWithLength<Derived>::operator=(const AbstractArrayWithLength& array) {
        if (this != &array) {
            length = array.length;
            Base::operator=(array);
        }
        return *this;
    }

    template<class Derived>
    AbstractArrayWithLength<Derived>&
    AbstractArrayWithLength<Derived>::operator=(AbstractArrayWithLength&& array) noexcept {
        length = array.length;
        Base::operator=(std::move(array));
        return *this;
    }
    /**
     * Get the last element in the array and remove it from the array.
     */
    template<class Derived>
    typename AbstractArrayWithLength<Derived>::T AbstractArrayWithLength<Derived>::cutLast() {
        assert(length > 0);
        --length;
        if(QTypeInfo<T>::isComplex)
            return T(std::move(arr[length]));
        else
            return arr[length];
    }
    /**
     * Low level api. Designed for performance.
     * Allocate a element at the end and increase the length.
     * This function can be used when you are sure the current capacity is enough.
     */
    template<class Derived>
    inline void AbstractArrayWithLength<Derived>::grow(const T& t) {
        assert(length < Base::getDerived().getCapacity());
        Base::allocate(t, length++);
    }

    template<class Derived>
    inline void AbstractArrayWithLength<Derived>::grow(T&& t) {
        assert(length < Base::getDerived().getCapacity());
        Base::allocate(std::move(t), length++);
    }

    template<class Derived>
    void AbstractArrayWithLength<Derived>::removeAt(size_t index) {
        assert(index < length);
        if(QTypeInfo<T>::isComplex)
            (arr + index)->~T();
        --length;
        memmove(arr + index, arr + index + 1, (length - index) * sizeof(T));
    }

    template<class Derived>
    void AbstractArrayWithLength<Derived>::clear() noexcept {
        for (size_t i = 0; i < length; ++i)
            (arr + i)->~T();
        length = 0;
    }

    template<class Derived>
    void AbstractArrayWithLength<Derived>::insert(const T& t, size_t index) {
        assert(length < Base::getCapacity());
        memmove(arr + index + 1, arr + index, length - index);
        Base::init(t, index);
        Base::setLength(length + 1);
    }

    template<class Derived>
    void AbstractArrayWithLength<Derived>::insert(T&& t, size_t index) {
        assert(length < Base::getCapacity());
        memmove(arr + index + 1, arr + index, length - index);
        Base::init(std::move(t), index);
        Base::setLength(length + 1);
    }

    template<class Derived>
    void AbstractArrayWithLength<Derived>::swap(AbstractArrayWithLength& array) {
        Base::swap(array);
        std::swap(length, array.length);
    }
}