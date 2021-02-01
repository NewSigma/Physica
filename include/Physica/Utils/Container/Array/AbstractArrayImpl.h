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

#include <cstring>

namespace Physica::Utils::Internal {
    template<class Pointer, class Derived>
    Iterator<Pointer, AbstractArray<Derived>>&
    Iterator<Pointer, AbstractArray<Derived>>::operator=(const Iterator& ite) { //NOLINT Self assign is ok.
        p = ite.p;
        return *this;
    }

    template<class Pointer, class Derived>
    Iterator<Pointer, AbstractArray<Derived>>&
    Iterator<Pointer, AbstractArray<Derived>>::operator++() {
        ++p;
        return *this;
    }

    template<class Pointer, class Derived>
    const Iterator<Pointer, AbstractArray<Derived>>
    Iterator<Pointer, AbstractArray<Derived>>::operator++(int) {
        return Iterator(p++);
    }

    template<class Pointer, class Derived>
    ReverseIterator<Pointer, AbstractArray<Derived>>&
    ReverseIterator<Pointer, AbstractArray<Derived>>::operator=(const ReverseIterator& ite) { //NOLINT Self assign is ok.
        p = ite.p;
        return *this;
    }

    template<class Pointer, class Derived>
    ReverseIterator<Pointer, AbstractArray<Derived>>&
    ReverseIterator<Pointer, AbstractArray<Derived>>::operator++() {
        --p;
        return *this;
    }

    template<class Pointer, class Derived>
    const ReverseIterator<Pointer, AbstractArray<Derived>>
    ReverseIterator<Pointer, AbstractArray<Derived>>::operator++(int) {
        return ReverseIterator(p--);
    }

    template<class Derived>
    AbstractArray<Derived>::AbstractArray(size_t capacity)
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))) {}

    template<class Derived>
    AbstractArray<Derived>::AbstractArray(const AbstractArray<Derived>& array)
            : AbstractArray(array.getDerived().getCapacity()) {
        const size_t length = array.getDerived().getLength();
        if(QTypeInfo<T>::isComplex)
            for(size_t i = 0; i < length; ++i)
                allocate(array[i], i);
        else
            memcpy(arr, array.arr, length * sizeof(T));
    }

    template<class Derived>
    AbstractArray<Derived>::AbstractArray(AbstractArray<Derived>&& array) noexcept : arr(array.arr) {
        array.arr = nullptr;
    }

    template<class Derived>
    AbstractArray<Derived>::~AbstractArray() {
        //Elements should call ~T() at subclasses.
        free(arr);
    }

    template<class Derived>
    AbstractArray<Derived>& AbstractArray<Derived>::operator=(const AbstractArray& array) {
        if (this != &array) {
            getDerived().~Derived();
            arr = reinterpret_cast<T*>(malloc(array.getDerived().getCapacity() * sizeof(T)));
            const size_t length = array.getDerived().getLength();
            if(QTypeInfo<T>::isComplex)
                for(size_t i = 0; i < length; ++i)
                    allocate(array[i], i);
            else
                memcpy(arr, array.arr, length * sizeof(T));
        }
        return *this;
    }

    template<class Derived>
    AbstractArray<Derived>& AbstractArray<Derived>::operator=(AbstractArray&& array) noexcept {
        getDerived().~Derived();
        arr = array.arr;
        array.arr = nullptr;
        return *this;
    }

    template<class Derived>
    inline typename AbstractArray<Derived>::T& AbstractArray<Derived>::operator[](size_t index) {
        Q_ASSERT(index < getDerived().getLength());
        return arr[index];
    }

    template<class Derived>
    inline const typename AbstractArray<Derived>::T& AbstractArray<Derived>::operator[](size_t index) const {
        Q_ASSERT(index < getDerived().getLength());
        return arr[index];
    }

    template<class Derived>
    bool AbstractArray<Derived>::operator==(const AbstractArray& array) const {
        if (getDerived().getLength() != array.getDerived().getLength())
            return false;
        if (getDerived().getCapacity() != array.getDerived().getLength())
            return false;
        for (size_t i = 0; i < getDerived().getLength(); ++i)
            if (operator[](i) != array[i])
                return false;
        return true;
    }
    /**
     * Low level api. Designed for performance.
     * Simply allocate a T at position \param index.
     * You must ensure position \index is usable or a memory leak will occur.
     * 
     * Expose this function or not depends on the subclasses.
     */
    template<class Derived>
    inline void AbstractArray<Derived>::allocate(const T& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(t);
        else
            *(arr + index) = t;
    }

    template<class Derived>
    inline void AbstractArray<Derived>::allocate(T&& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(std::move(t));
        else
            *(arr + index) = t;
    }

    template<class Derived>
    inline void AbstractArray<Derived>::swap(AbstractArray<Derived>& array) noexcept {
        std::swap(arr, array.arr);
    }
}
