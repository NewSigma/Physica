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

namespace Physica::Utils {
    template<class T> class DeviceAllocator;
}

namespace Physica::Utils::Internal {
    template<class Pointer, class Derived, class Allocator>
    __host__ __device__ ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>&
    ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>::operator=(const ContainerIterator& ite) { //NOLINT Self assign is ok.
        p = ite.p;
        return *this;
    }

    template<class Pointer, class Derived, class Allocator>
    __host__ __device__ ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>
    ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>::operator+(difference_type n) const {
        return ContainerIterator(p + n);
    }

    template<class Pointer, class Derived, class Allocator>
    __host__ __device__ ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>
    ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>::operator-(difference_type n) const {
        return ContainerIterator(p - n);
    }

    template<class Pointer, class Derived, class Allocator>
    __host__ __device__ ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>&
    ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>::operator++() {
        ++p;
        return *this;
    }

    template<class Pointer, class Derived, class Allocator>
    __host__ __device__ const ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>
    ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>::operator++(int) {
        return ContainerIterator(p++);
    }

    template<class Pointer, class Derived, class Allocator>
    __host__ __device__ ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>&
    ContainerIterator<Pointer, ArrayBase<Derived, Allocator>>::operator--() {
        --p;
        return *this;
    }

    template<class Pointer, class Derived, class Allocator>
    __host__ __device__ ReverseContainerIterator<Pointer, ArrayBase<Derived, Allocator>>&
    ReverseContainerIterator<Pointer, ArrayBase<Derived, Allocator>>::operator=(const ReverseContainerIterator& ite) { //NOLINT Self assign is ok.
        p = ite.p;
        return *this;
    }

    template<class Pointer, class Derived, class Allocator>
    __host__ __device__ ReverseContainerIterator<Pointer, ArrayBase<Derived, Allocator>>&
    ReverseContainerIterator<Pointer, ArrayBase<Derived, Allocator>>::operator++() {
        --p;
        return *this;
    }

    template<class Pointer, class Derived, class Allocator>
    __host__ __device__ const ReverseContainerIterator<Pointer, ArrayBase<Derived, Allocator>>
    ReverseContainerIterator<Pointer, ArrayBase<Derived, Allocator>>::operator++(int) {
        return ReverseContainerIterator(p--);
    }

    template<class Derived, class Allocator>
    __host__ __device__ inline typename ArrayBase<Derived, Allocator>::LValueReferenceType
    ArrayBase<Derived, Allocator>::operator[](size_t index) {
        assert(index < Base::getDerived().getLength());
        return data()[index];
    }

    template<class Derived, class Allocator>
    __host__ __device__ inline typename ArrayBase<Derived, Allocator>::ConstLValueReferenceType
    ArrayBase<Derived, Allocator>::operator[](size_t index) const {
        assert(index < Base::getDerived().getLength());
        return data()[index];
    }

    template<class Derived, class Allocator>
    bool ArrayBase<Derived, Allocator>::operator==(const ArrayBase& array) const {
        if (Base::getDerived().getLength() != array.getDerived().getLength())
            return false;
        if (Base::getDerived().getCapacity() != array.getDerived().getLength())
            return false;
        for (size_t i = 0; i < Base::getDerived().getLength(); ++i)
            if (operator[](i) != array[i])
                return false;
        return true;
    }
}
