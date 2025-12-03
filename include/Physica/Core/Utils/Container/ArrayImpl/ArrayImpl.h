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
    __host__ __device__ Array<T, Length, Allocator>::Array([[maybe_unused]] size_t length, auto&&... args) {
        assert(length == Length);
        if constexpr (!Base::template isTrivialDefaultConstruct<decltype(args)...>()) {
            for (size_t i = 0; i < Length; ++i)
                arr[i] = T(std::forward<decltype(args)>(args)...);
        }
    }

    template<class T, size_t Length, class Allocator>
    __host__ __device__ Array<T, Length, Allocator>::Array(std::initializer_list<T> list) {
        assert(list.size() == Length);
        unsigned int i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            arr[i] = *ite;
    }

    template<class T, size_t Length, class Allocator>
    template<size_t OtherLength, class OtherAlloc>
    Array<T, Length, Allocator>::Array(const Array<T, OtherLength, OtherAlloc>& other) noexcept : Array(read(Length, other.data())) {
        static_assert(OtherLength == Length || OtherLength == Dynamic, "[Error]: Length do not match");
        assert(other.getLength() >= Length);
    }
    /**
     * Initializing new elements will not work. A fixed array is assumed to be initialized upon construction.
     */
    template<class T, size_t Length, class Allocator>
    __host__ __device__ void Array<T, Length, Allocator>::resize([[maybe_unused]] size_t length, auto&&...) noexcept {
        assert(length == Length);
    }

    template<class T, size_t Length, class Allocator>
    __host__ __device__ void Array<T, Length, Allocator>::zeros() noexcept {
        static_assert(std::is_trivially_copyable<T>::value, "[Error]: zeros() does not apply to non-trivial type");
        memset(arr.data(), 0, Length * sizeof(T));
    }

    template<class T, size_t Length, class Allocator>
    __host__ __device__ void Array<T, Length, Allocator>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        #ifdef __CUDA_ARCH__
            thrust::swap(arr, obj.arr);
        #else
            arr.swap(obj.arr);
        #endif
    }
    /**
     * Helper function that communicates with C libraries.
     */
    template<class T, size_t Length, class Allocator>
    __host__ __device__ auto Array<T, Length, Allocator>::read([[maybe_unused]] size_t length, const T* __restrict p) noexcept -> This {
        static_assert(std::is_trivially_copyable<T>::value, "[Error]: C type must be trivial");
        assert(length == Length && "[Error]: Length do not match");
        assert(p != nullptr);
        This result(Length);
        std::memcpy(result.arr.data(), p, Length * sizeof(ElemType));
        return result;
    }

    template<class T, size_t Length, class Allocator>
    size_t Array<T, Length, Allocator>::toIndex1D(const IndexType& __restrict shape, const IndexType& __restrict indices) noexcept {
        size_t index = 0;
        size_t stride = 1;
        for (int i = static_cast<int>(shape.getLength()) - 1; i >= 0; --i) {
            assert(indices[i] < shape[i] && "[Error]: Index out of range");
            index += indices[i] * stride;
            stride *= shape[i];
        }
        return index;
    }

    template<class T, size_t Length, class Allocator>
    auto Array<T, Length, Allocator>::toIndexND(const IndexType& shape, size_t index) noexcept -> IndexType {
        const int dim = shape.getLength();
        IndexType indices(shape.size());
        size_t remaining = index;
        for (int i = 0; i < dim; ++i) {
            size_t stride = 1;
            for (int j = i + 1; j < dim; ++j)
                stride *= shape[j]; 
            indices[i] = remaining / stride;
            assert(indices[i] < shape[i] && "[Error]: Index out of range");
            remaining %= stride;
        }
        return indices;
    }
    ///////////////////////////////////////Array<T, Dynamic, Allocator>//////////////////////////////////////////
    template<class T, class Allocator>
    Array<T, Dynamic, Allocator>::Array(size_t length_, auto&&... args) noexcept(std::is_nothrow_constructible<T, decltype(args)...>::value)
            : length(length_), capacity(length_), alloc() {
        if (capacity == 0)
            return;

        arr = alloc.allocate(capacity);
        if constexpr (!Base::template isTrivialDefaultConstruct<decltype(args)...>()) {
            for (size_t i = 0; i < length_; ++i)
                alloc.construct(arr + i, std::forward<decltype(args)>(args)...);
        }
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Allocator>::Array(std::initializer_list<T> list) noexcept : Array(list.size()) {
        size_t i = 0;
        const auto end = list.end();
        for (auto ite = list.begin(); ite != end; ++ite, ++i)
            alloc.construct(arr + i, std::move(*ite));
    }

    template<class T, class Allocator>
    template<size_t Length, class OtherAlloc>
    Array<T, Dynamic, Allocator>::Array(const Array<T, Length, OtherAlloc>& other) noexcept(std::is_nothrow_copy_assignable<T>::value)
            : Array(other.getLength()) {
        for (size_t i = 0; i < getLength(); ++i)
            arr[i] = other[i];
    }

    template<class T, class Allocator>
    Array<T, Dynamic, Allocator>::Array(const This& obj) noexcept(std::is_nothrow_copy_constructible<T>::value)
            : length(obj.length), capacity(obj.capacity), alloc() {
        if (capacity == 0)
            return;

        arr = alloc.allocate(capacity);
        if constexpr (!std::is_trivially_copyable<ElemType>::value)
            for(size_t i = 0; i < length; ++i)
                alloc.construct(arr + i, obj[i]);
        else
            memcpy(arr, obj.arr, length * sizeof(ElemType));
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
        if constexpr (!std::is_trivially_copyable<ElemType>::value)
            if (arr != nullptr)
                for(size_t i = 0; i < length; ++i)
                    alloc.destroy(arr + i);
        alloc.deallocate(arr, length);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::grow(auto&&... args) noexcept {
        assert(length < getCapacity() && "[Error]: You must make sure capacity is enough before calling grow()");
        alloc.construct(arr + length++, std::forward<decltype(args)>(args)...);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::append(auto&&... args) noexcept {
        if (length == capacity) [[unlikely]]
            doubleSpace();
        grow(std::forward<decltype(args)>(args)...);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::insert(size_t index, auto&&... args) noexcept {
        assert(index <= length);
        if (length == capacity)
            doubleSpace();
        memmove(arr + index + 1, arr + index, (length - index) * sizeof(T));
        alloc.construct(arr + index, std::forward<decltype(args)>(args)...);
        setLength(length + 1);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::reserve(size_t size) noexcept {
        if (size > getCapacity()) {
            arr = alloc.reallocate(arr, size, capacity);
            capacity = size;
        }
    }
    /**
     * resize() is inlined, if there is no need to resize, we avoid the overhead of slow resizeImpl() call
     */
    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::resize(size_t size, auto&&... args) noexcept {
        if (length != size)
            resizeImpl(size, std::forward<decltype(args)>(args)...);
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::squeeze() noexcept {
        if (length == 0)
            *this = Array();
        else {
            arr = alloc.reallocate(arr, length, capacity);
            capacity = length;
        }
    }
    /*!
     * Increase the capacity.
     * This function can be used when you are sure the new \param size is larger than the old capacity.
     */
    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::increase(size_t size) noexcept {
        assert(size >= capacity);
        arr = alloc.reallocate(arr, size, capacity);
        capacity = size;
    }
    /*!
     * Decrease the capacity.
     * This function can be used when you are sure the new \param size is shorter than the old capacity.
     */
    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::decrease(size_t size) noexcept {
        assert(size <= capacity);
        if(!std::is_trivially_copyable<T>::value) {
            for(size_t i = size; i < length; ++i)
                (arr + i)->~T();
        }
        arr = alloc.reallocate(arr, size);
        length = capacity = size;
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::clear() noexcept {
        for (size_t i = 0; i < length; ++i)
            alloc.destroy(arr + i);
        length = 0;
    }

    template<class T, class Allocator>
    auto Array<T, Dynamic, Allocator>::release() noexcept -> pointer {
        pointer p = arr;
        arr = nullptr;
        length = capacity = 0;
        return p;
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::doubleSpace() noexcept {
        increase((capacity * 2) + ((MinDeltaSpace + sizeof(T) - 1) / sizeof(T)));
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::zeros() noexcept {
        static_assert(std::is_trivially_copyable<T>::value, "[Error]: zeros() does not apply to non-trivial type");
        memset(arr, 0, length * sizeof(T));
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        // FIXME: Make use of trivially_relocatable once we dump to CXX26
        std::swap(arr, obj.arr);
        std::swap(length, obj.length);
        std::swap(capacity, obj.capacity);
    }

    template<class T, class Allocator>
    __host__ __device__ auto Array<T, Dynamic, Allocator>::data() noexcept -> pointer {
        if constexpr (Align == Dynamic)
            return arr;
        else
            return std::assume_aligned<Align, ElemType>(arr);
    }

    template<class T, class Allocator>
    __host__ __device__ auto Array<T, Dynamic, Allocator>::data() const noexcept -> const_pointer {
        return const_cast<This&>(*this).data();
    }
    /**
     * Low level api. Designed for performance. Elements between old length and \param size have not allocated. DO NOT try to visit them.
     */
    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::setLength(size_t size) noexcept {
        assert(size <= getCapacity() && "[Error]: Requiring more elements than the array have");
        if constexpr (!std::is_trivially_copyable<ElemType>::value) {
            assert(length <= size && "[Error]: setLength() cannot destruct unused elements, memory leak is expected");
        }
        length = size;
    }
    /**
     * Helper function that communicates with C libraries.
     */
    template<class T, class Allocator>
    auto Array<T, Dynamic, Allocator>::read(size_t length, const T* __restrict p) -> This {
        static_assert(std::is_trivially_copyable<T>::value, "[Error]: C type must be trivial");
        assert(p != nullptr);
        This result(length);
        memcpy(result.arr, p, length * sizeof(ElemType));
        return result;
    }

    template<class T, class Allocator>
    void Array<T, Dynamic, Allocator>::resizeImpl(size_t size, auto&&... args) noexcept {
        if (capacity < size)
            reserve(size);

        if (length > size) {
            if constexpr (!std::is_trivially_copyable<T>::value)
                for (size_t i = size; i < length; ++i)
                    (arr + i)->~T();
            length = size;
        }
        else {
            if constexpr (Base::template isTrivialDefaultConstruct<decltype(args)...>())
                length = size;
            else {
                for (; length < size; ++length)
                    alloc.construct(arr + length, std::forward<decltype(args)>(args)...);
            }
        }
    }
}
