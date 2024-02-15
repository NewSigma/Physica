/*
 * Copyright 2020-2024 WeiBo He.
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
    template<class Derived, class Allocator>
    DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(size_t capacity)
            : alloc(), length(0) {
        arr = alloc.allocate(capacity);
    }

    template<class Derived, class Allocator>
    DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(size_t length_, size_t capacity)
            : DynamicArrayBase(capacity) {
        assert(length <= capacity);
        length = length_;
    }

    template<class Derived, class Allocator>
    DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(size_t length_, pointer arr_)
            : arr(arr_), alloc(), length(length_) {}

    template<class Derived, class Allocator>
    DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(
            const DynamicArrayBase& array) : DynamicArrayBase(array.getDerived().getCapacity()) {
        length = array.length;
        if constexpr (!std::is_trivial<ValueType>::value)
            for(size_t i = 0; i < length; ++i)
                alloc.construct(arr + i, array[i]);
        else
            memcpy(arr, array.arr, length * sizeof(ValueType));
    }
    
    template<class Derived, class Allocator>
    DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(
            DynamicArrayBase&& array) noexcept : arr(array.arr), alloc(), length(array.length) {
        array.arr = nullptr;
        array.length = 0;
    }
    
    template<class Derived, class Allocator>
    DynamicArrayBase<Derived, Allocator>::~DynamicArrayBase() {
        if constexpr (!std::is_trivial<ValueType>::value)
            if (arr != nullptr)
                for(size_t i = 0; i < length; ++i)
                    alloc.destroy(arr + i);
        alloc.deallocate(arr, length);
    }

    template<class Derived, class Allocator>
    DynamicArrayBase<Derived, Allocator>& DynamicArrayBase<Derived, Allocator>::operator=(DynamicArrayBase array) noexcept {
        swap(array);
        return *this;
    }
    /**
     * Get the last element in the array and remove it from the array.
     */
    template<class Derived, class Allocator>
    typename DynamicArrayBase<Derived, Allocator>::ValueType DynamicArrayBase<Derived, Allocator>::cutLast() {
        assert(length > 0);
        --length;
        if constexpr (!std::is_trivial<ValueType>::value)
            return T(std::move(arr[length]));
        else
            return arr[length];
    }

    template<class Derived, class Allocator>
    template<class... Args>
    inline void DynamicArrayBase<Derived, Allocator>::grow(Args&&... args) {
        assert(length < Base::getDerived().getCapacity() && "[Error]: You must make sure capacity is enough before calling grow()");
        alloc.construct(arr + length++, std::forward<Args>(args)...);
    }

    template<class Derived, class Allocator>
    void DynamicArrayBase<Derived, Allocator>::removeAt(size_t index) {
        assert(index < length);
        if constexpr (!std::is_trivial<ValueType>::value)
            alloc.destroy(arr + index);
        --length;
        memmove(arr + index, arr + index + 1, (length - index) * sizeof(ValueType));
    }

    template<class Derived, class Allocator>
    void DynamicArrayBase<Derived, Allocator>::clear() noexcept {
        for (size_t i = 0; i < length; ++i)
            alloc.destroy(arr + i);
        length = 0;
    }
    /**
     * Low level api. Designed for performance. Elements between old length and \param size have not allocated. DO NOT try to visit them.
     */
    template<class Derived, class Allocator>
    inline void DynamicArrayBase<Derived, Allocator>::setLength(size_t size) {
        assert(size <= Base::getDerived().getCapacity() && "[Error]: Requiring more elements than the array have");
        if constexpr (!std::is_trivial<ValueType>::value) {
            assert(length <= size && "[Error]: setLength() cannot destruct unused elements, memory leak is expected");
        }
        length = size;
    }

    template<class Derived, class Allocator>
    void DynamicArrayBase<Derived, Allocator>::swap(DynamicArrayBase& __restrict array) noexcept {
        assert(this != &array && "[Error]: Self swap is likely a bug");
        std::swap(arr, array.arr);
        std::swap(length, array.length);
    }

    template<class Derived, class Allocator>
    inline typename DynamicArrayBase<Derived, Allocator>::pointer
    DynamicArrayBase<Derived, Allocator>::release() noexcept {
        pointer copy = arr;
        length = 0;
        arr = nullptr;
        return copy;
    }
}