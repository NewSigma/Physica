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

namespace Physica::Utils::Internal {
    template<class Derived>
    DynamicArrayBase<Derived>::DynamicArrayBase(size_t capacity)
        : Base(capacity), length(0) {}

    template<class Derived>
    DynamicArrayBase<Derived>::DynamicArrayBase(size_t length_, size_t capacity)
            : Base(capacity), length(length_) {
        assert(length <= capacity);
    }

    template<class Derived>
    DynamicArrayBase<Derived>::DynamicArrayBase(size_t length_, T* __restrict arr_)
            : Base(arr_), length(length_) {}

    template<class Derived>
    DynamicArrayBase<Derived>::DynamicArrayBase(
        const DynamicArrayBase& array) : Base(array), length(array.length) {}
    
    template<class Derived>
    DynamicArrayBase<Derived>::DynamicArrayBase(
        DynamicArrayBase&& array) noexcept : Base(std::move(array)), length(array.length) {}
    /**
     * Get the last element in the array and remove it from the array.
     */
    template<class Derived>
    typename DynamicArrayBase<Derived>::T DynamicArrayBase<Derived>::cutLast() {
        assert(length > 0);
        --length;
        if constexpr (!std::is_trivial<T>::value)
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
    inline void DynamicArrayBase<Derived>::grow(const T& t) {
        assert(length < Base::getDerived().getCapacity());
        alloc.construct(arr + length++, t);
    }

    template<class Derived>
    inline void DynamicArrayBase<Derived>::grow(T&& t) {
        assert(length < Base::getDerived().getCapacity());
        alloc.construct(arr + length++, std::move(t));
    }

    template<class Derived>
    void DynamicArrayBase<Derived>::removeAt(size_t index) {
        assert(index < length);
        if constexpr (!std::is_trivial<T>::value)
            (arr + index)->~T();
        --length;
        memmove(arr + index, arr + index + 1, (length - index) * sizeof(T));
    }

    template<class Derived>
    void DynamicArrayBase<Derived>::clear() noexcept {
        for (size_t i = 0; i < length; ++i)
            (arr + i)->~T();
        length = 0;
    }

    template<class Derived>
    void DynamicArrayBase<Derived>::insert(const T& t, size_t index) {
        assert(length < Base::getCapacity());
        memmove(arr + index + 1, arr + index, length - index);
        alloc.construct(arr + index, t);
        Base::setLength(length + 1);
    }

    template<class Derived>
    void DynamicArrayBase<Derived>::insert(T&& t, size_t index) {
        assert(length < Base::getCapacity());
        memmove(arr + index + 1, arr + index, length - index);
        alloc.construct(arr + index, std::move(t));
        Base::setLength(length + 1);
    }

    template<class Derived>
    void DynamicArrayBase<Derived>::swap(DynamicArrayBase& array) {
        Base::swap(array);
        std::swap(length, array.length);
    }
}