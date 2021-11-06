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
    template<class Derived, class Allocator>
    __host__ __device__ DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(size_t capacity)
            : alloc(), length(0) {
        arr = alloc.allocate(capacity);
    }

    template<class Derived, class Allocator>
    __host__ __device__ DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(size_t length_, size_t capacity)
            : DynamicArrayBase(capacity) {
        assert(length <= capacity);
        length = length_;
    }

    template<class Derived, class Allocator>
    __host__ __device__ DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(size_t length_, PointerType arr_)
            : arr(arr_), alloc(), length(length_) {}

    template<class Derived, class Allocator>
    __host__ __device__ DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(
            const DynamicArrayBase& array) : DynamicArrayBase(array.getDerived().getCapacity()) {
        length = array.length;
        if constexpr (!std::is_trivial<ValueType>::value)
            for(size_t i = 0; i < length; ++i)
                AllocatorTraits::construct(alloc, arr + i, array[i]);
        else {
        #ifdef PHYSICA_CUDA
            #ifdef __CUDA__ARCH__
                memcpy(arr, array.arr, length * sizeof(ValueType));
            #else
                if constexpr (std::is_same<allocator_type, DeviceAllocator<ValueType>>::value)
                    cudaMemcpy(arr.get(), array.arr.get(), length * sizeof(ValueType), cudaMemcpyDeviceToDevice);
                else
                    memcpy(arr, array.arr, length * sizeof(ValueType));
            #endif
        #else
            memcpy(arr, array.arr, length * sizeof(ValueType));
        #endif
        }
    }
    
    template<class Derived, class Allocator>
    __host__ __device__ DynamicArrayBase<Derived, Allocator>::DynamicArrayBase(
            DynamicArrayBase&& array) noexcept : arr(array.arr), alloc(), length(array.length) {
        array.arr = nullptr;
        array.length = 0;
    }
    
    template<class Derived, class Allocator>
    __host__ __device__ DynamicArrayBase<Derived, Allocator>::~DynamicArrayBase() {
        if constexpr (!std::is_trivial<ValueType>::value)
            if (arr != nullptr)
                for(size_t i = 0; i < length; ++i)
                    AllocatorTraits::destroy(alloc, arr + i);
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
    /**
     * Low level api. Designed for performance.
     * Allocate a element at the end and increase the length.
     * This function can be used when you are sure the current capacity is enough.
     */
    template<class Derived, class Allocator>
    __host__ __device__ inline void DynamicArrayBase<Derived, Allocator>::grow(ConstLValueReferenceType t) {
        assert(length < Base::getDerived().getCapacity());
        alloc.construct(arr + length++, t);
    }

    template<class Derived, class Allocator>
    __host__ __device__ inline void DynamicArrayBase<Derived, Allocator>::grow(RValueReferenceType t) {
        assert(length < Base::getDerived().getCapacity());
        alloc.construct(arr + length++, std::move(t));
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
    __host__ __device__ void DynamicArrayBase<Derived, Allocator>::clear() noexcept {
        for (size_t i = 0; i < length; ++i)
            alloc.destroy(arr + i);
        length = 0;
    }

    template<class Derived, class Allocator>
    void DynamicArrayBase<Derived, Allocator>::insert(ConstLValueReferenceType t, size_t index) {
        assert(length < Base::getCapacity());
        memmove(arr + index + 1, arr + index, length - index);
        alloc.construct(arr + index, t);
        Base::setLength(length + 1);
    }

    template<class Derived, class Allocator>
    void DynamicArrayBase<Derived, Allocator>::insert(RValueReferenceType t, size_t index) {
        assert(length < Base::getCapacity());
        memmove(arr + index + 1, arr + index, length - index);
        alloc.construct(arr + index, std::move(t));
        Base::setLength(length + 1);
    }

    template<class Derived, class Allocator>
    __host__ __device__ void DynamicArrayBase<Derived, Allocator>::swap(DynamicArrayBase& array) {
        std::swap(arr, array.arr);
        std::swap(length, array.length);
    }
}