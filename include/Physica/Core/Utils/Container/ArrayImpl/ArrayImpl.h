/*
 * Copyright 2020-2025 Weibo He.
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
#ifdef PHYSICA_CUDA
    #include <thrust/swap.h>
#endif
#include "../Array.h"

namespace Physica {
    //////////////////////////////////////////Array<T, Length, Allocator>//////////////////////////////////////////
    template<class T, size_t Length, class Allocator>
    template<class... Args>
    __host__ __device__ Array<T, Length, Allocator>::Array([[maybe_unused]] size_t length, Args&&... args) {
        assert(length == Length);
        if constexpr (!Base::template isTrivialDefaultConstruct<Args...>()) {
            for (size_t i = 0; i < Length; ++i)
                *(arr + i) = T(std::forward<Args>(args)...);
        }
    }

    template<class T, size_t Length, class Allocator>
    __host__ __device__ Array<T, Length, Allocator>::Array(std::initializer_list<T> list) {
        assert(list.size() <= Length);
        unsigned int i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            *(arr + i) = *ite;
    }

    template<class T, size_t Length, class Allocator>
    template<class OtherAlloc>
    Array<T, Length, Allocator>::Array(const Array<T, Length, OtherAlloc>& other) {
        for (size_t i = 0; i < Length; ++i)
            arr[i] = other[i];
    }
    /**
     * Initializing new elements will not work. A fixed array is assumed to be initialized upon construction.
     */
    template<class T, size_t Length, class Allocator>
    template<class... Args>
    __host__ __device__ inline void Array<T, Length, Allocator>::resize([[maybe_unused]] size_t length, Args&&...) {
        assert(length == Length);
    }

    template<class T, size_t Length, class Allocator>
    __host__ __device__ void Array<T, Length, Allocator>::swap(Array& __restrict array) noexcept {
        assert(this != &array && "[Error]: Self swap is likely a bug");
        for (size_t i = 0; i < Length; ++i) {
        #ifdef __CUDA_ARCH__
            thrust::swap(arr[i], array[i]);
        #else
            std::swap(arr[i], array[i]);
        #endif
        }
    }
    ///////////////////////////////////////Array<T, Dynamic, Allocator>//////////////////////////////////////////
    template<class T, class Allocator>
    Array<T, Dynamic, Allocator>::Array() : arr(nullptr), length(0), capacity(0), alloc() {}

    template<class T, class Allocator>
    template<class... Args>
    Array<T, Dynamic, Allocator>::Array(size_t length_, Args&&... args) : length(length_), capacity(length_), alloc() {
        arr = alloc.allocate(capacity);
        if constexpr (!Base::template isTrivialDefaultConstruct<Args...>()) {
            for (size_t i = 0; i < length_; ++i)
                alloc.construct(arr + i, std::forward<Args>(args)...);
        }
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Allocator>::Array(std::initializer_list<T> list) : Array(list.size()) {
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            alloc.construct(arr + i, *ite);
        setLength(list.size());
    }

    template<class T, class Allocator>
    template<class OtherAlloc>
    Array<T, Dynamic, Allocator>::Array(const Array<T, Dynamic, OtherAlloc>& other) : Array(other.getLength()) {
        for (size_t i = 0; i < getLength(); ++i)
            arr[i] = other[i];
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Allocator>::Array(const This& array) : length(array.length), capacity(array.capacity), alloc() {
        arr = alloc.allocate(capacity);
        if constexpr (!std::is_trivial<ElemType>::value)
            for(size_t i = 0; i < length; ++i)
                alloc.construct(arr + i, array[i]);
        else
            memcpy(arr, array.arr, length * sizeof(ElemType));
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Allocator>::Array(This&& array) noexcept
            : arr(array.arr), length(array.length), capacity(array.capacity), alloc() {
        array.arr = nullptr;
        array.length = 0;
        array.capacity = 0;
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Allocator>::~Array() {
        if constexpr (!std::is_trivial<ElemType>::value)
            if (arr != nullptr)
                for(size_t i = 0; i < length; ++i)
                    alloc.destroy(arr + i);
        alloc.deallocate(arr, length);
    }

    template<class T, class Allocator>
    template<class... Args>
    inline void Array<T, Dynamic, Allocator>::grow(Args&&... args) {
        assert(length < getCapacity() && "[Error]: You must make sure capacity is enough before calling grow()");
        alloc.construct(arr + length++, std::forward<Args>(args)...);
    }

    template<class T, class Allocator>
    template<class... Args>
    inline void Array<T, Dynamic, Allocator>::append(Args&&... args) {
        if (length == capacity) [[unlikely]]
            doubleSpace();
        grow(std::forward<Args>(args)...);
    }

    template<class T, class Allocator>
    template<class... Args>
    void Array<T, Dynamic, Allocator>::insert(size_t index, Args&&... args) {
        assert(index <= length);
        if (length == capacity)
            doubleSpace();
        memmove(arr + index + 1, arr + index, (length - index) * sizeof(T));
        alloc.construct(arr + index, std::forward<Args>(args)...);
        setLength(length + 1);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::reserve(size_t size) {
        if (size > getCapacity()) {
            arr = alloc.reallocate(arr, size, capacity);
            capacity = size;
        }
    }
    /**
     * resize() is inlined, if there is no need to resize, we avoid the overhead of a slow function call
     */
    template<class T, class Allocator>
    template<class... Args>
    inline void Array<T, Dynamic, Allocator>::resize(size_t size, Args&&... args) {
        if (length != size)
            resizeImpl<Args...>(size, std::forward<Args>(args)...);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::squeeze() {
        arr = alloc.reallocate(arr, length, capacity);
        capacity = length;
    }
    /*!
     * Increase the capacity.
     * This function can be used when you are sure the new \param size is larger than the old capacity.
     */
    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::increase(size_t size) {
        assert(size >= capacity);
        arr = alloc.reallocate(arr, size, capacity);
        capacity = size;
    }
    /*!
     * Decrease the capacity.
     * This function can be used when you are sure the new \param size is shorter than the old capacity.
     */
    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::decrease(size_t size) {
        assert(size <= capacity);
        if(!std::is_trivial<T>::value) {
            for(size_t i = size; i < length; ++i)
                (arr + i)->~T();
        }
        arr = alloc.reallocate(arr, size);
        length = capacity = size;
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::swap(Array& __restrict array) noexcept {
        assert(this != &array && "[Error]: Self swap is likely a bug");
        std::swap(arr, array.arr);
        std::swap(length, array.length);
        std::swap(capacity, array.capacity);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::clear() noexcept {
        for (size_t i = 0; i < length; ++i)
            alloc.destroy(arr + i);
        length = 0;
    }

    template<class T, class Allocator>
    inline Array<T, Dynamic, Allocator>::pointer Array<T, Dynamic, Allocator>::release() noexcept {
        pointer p = arr;
        arr = nullptr;
        length = capacity = 0;
        return p;
    }
    /**
     * Low level api. Designed for performance. Elements between old length and \param size have not allocated. DO NOT try to visit them.
     */
    template<class T, class Allocator>
    inline void Array<T, Dynamic, Allocator>::setLength(size_t size) {
        assert(size <= getCapacity() && "[Error]: Requiring more elements than the array have");
        if constexpr (!std::is_trivial<ElemType>::value) {
            assert(length <= size && "[Error]: setLength() cannot destruct unused elements, memory leak is expected");
        }
        length = size;
    }

    template<class T, class Allocator>
    template<class... Args>
    void Array<T, Dynamic, Allocator>::resizeImpl(size_t size, Args&&... args) {
        if (capacity < size)
            reserve(size);

        if (length > size) {
            if constexpr (!std::is_trivial<T>::value)
                for (size_t i = size; i < length; ++i)
                    (arr + i)->~T();
            length = size;
        }
        else {
            if constexpr (Base::template isTrivialDefaultConstruct<Args...>())
                length = size;
            else {
                for (; length < size; ++length)
                    alloc.construct(arr + length, std::forward<Args>(args)...);
            }
        }
    }
}
