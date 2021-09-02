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
    __host__ __device__ ContainerIterator<Pointer, ArrayStorage<Derived>>&
    ContainerIterator<Pointer, ArrayStorage<Derived>>::operator=(const ContainerIterator& ite) { //NOLINT Self assign is ok.
        p = ite.p;
        return *this;
    }

    template<class Pointer, class Derived>
    __host__ __device__ ContainerIterator<Pointer, ArrayStorage<Derived>>&
    ContainerIterator<Pointer, ArrayStorage<Derived>>::operator++() {
        ++p;
        return *this;
    }

    template<class Pointer, class Derived>
    __host__ __device__ const ContainerIterator<Pointer, ArrayStorage<Derived>>
    ContainerIterator<Pointer, ArrayStorage<Derived>>::operator++(int) {
        return ContainerIterator(p++);
    }

    template<class Pointer, class Derived>
    __host__ __device__ ReverseContainerIterator<Pointer, ArrayStorage<Derived>>&
    ReverseContainerIterator<Pointer, ArrayStorage<Derived>>::operator=(const ReverseContainerIterator& ite) { //NOLINT Self assign is ok.
        p = ite.p;
        return *this;
    }

    template<class Pointer, class Derived>
    __host__ __device__ ReverseContainerIterator<Pointer, ArrayStorage<Derived>>&
    ReverseContainerIterator<Pointer, ArrayStorage<Derived>>::operator++() {
        --p;
        return *this;
    }

    template<class Pointer, class Derived>
    __host__ __device__ const ReverseContainerIterator<Pointer, ArrayStorage<Derived>>
    ReverseContainerIterator<Pointer, ArrayStorage<Derived>>::operator++(int) {
        return ReverseContainerIterator(p--);
    }

    template<class Derived>
    __host__ __device__ ArrayStorage<Derived>::ArrayStorage(size_t capacity)
            : alloc() {
        arr = alloc.allocate(capacity);
    }

    template<class Derived>
    __host__ __device__ ArrayStorage<Derived>::ArrayStorage(PointerType arr_) : arr(arr_), alloc() {}

    template<class Derived>
    __host__ __device__ ArrayStorage<Derived>::ArrayStorage(const ArrayStorage<Derived>& array)
            : ArrayStorage(array.getDerived().getCapacity()) {
        if constexpr (!std::is_trivial<ValueType>::value)
            for(size_t i = 0; i < array.getDerived().getLength(); ++i)
                AllocatorTraits::construct(alloc, arr + i, array[i]);
        else {
        #ifdef __CUDA__ARCH__
            memcpy(arr, array.arr, array.getDerived().getLength() * sizeof(ValueType));
        #else
            if constexpr (std::is_same<allocator_type, DeviceAllocator<ValueType>>::value)
                cudaMemcpy(arr.get(), array.arr.get(), array.getDerived().getLength() * sizeof(ValueType), cudaMemcpyDeviceToDevice);
            else
                memcpy(arr, array.arr, array.getDerived().getLength() * sizeof(ValueType));
        #endif
        }
    }

    template<class Derived>
    __host__ __device__ ArrayStorage<Derived>::ArrayStorage(ArrayStorage<Derived>&& array) noexcept
            : arr(array.arr)
            , alloc() {
        array.arr = nullptr;
    }

    template<class Derived>
    __host__ __device__ ArrayStorage<Derived>::~ArrayStorage() {
        const size_t length = Base::getDerived().getLength();
        if constexpr (!std::is_trivial<ValueType>::value)
            if (arr != nullptr)
                for(size_t i = 0; i < length; ++i)
                    AllocatorTraits::destroy(alloc, arr + i);
        alloc.deallocate(arr, length);
    }

    template<class Derived>
    __host__ __device__ inline typename ArrayStorage<Derived>::LValueReferenceType
    ArrayStorage<Derived>::operator[](size_t index) {
        assert(index < Base::getDerived().getLength());
        return arr[index];
    }

    template<class Derived>
    __host__ __device__ inline typename ArrayStorage<Derived>::ConstLValueReferenceType
    ArrayStorage<Derived>::operator[](size_t index) const {
        assert(index < Base::getDerived().getLength());
        return arr[index];
    }

    template<class Derived>
    bool ArrayStorage<Derived>::operator==(const ArrayStorage& array) const {
        if (Base::getDerived().getLength() != array.getDerived().getLength())
            return false;
        if (Base::getDerived().getCapacity() != array.getDerived().getLength())
            return false;
        for (size_t i = 0; i < Base::getDerived().getLength(); ++i)
            if (operator[](i) != array[i])
                return false;
        return true;
    }

    template<class Derived>
    __host__ __device__ inline void ArrayStorage<Derived>::swap(ArrayStorage<Derived>& array) noexcept {
        std::swap(arr, array.arr);
    }
}
