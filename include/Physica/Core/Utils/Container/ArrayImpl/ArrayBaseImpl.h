/*
 * Copyright 2020-2024 Weibo He.
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

namespace Physica::Core {
    template<class Derived, class Allocator>
    __host__ __device__ inline ArrayBase<Derived, Allocator>::lvalue_reference
    ArrayBase<Derived, Allocator>::operator[](size_t index) {
        assert(index < getLength() && "[Error]: Index overflow");
        return data()[index];
    }

    template<class Derived, class Allocator>
    __host__ __device__ inline ArrayBase<Derived, Allocator>::const_lvalue_reference
    ArrayBase<Derived, Allocator>::operator[](size_t index) const {
        assert(index < getLength() && "[Error]: Index overflow");
        return data()[index];
    }

    template<class Derived, class Allocator>
    bool ArrayBase<Derived, Allocator>::operator==(const ArrayBase& array) const {
        if (getLength() != array.getLength())
            return false;
        if (getCapacity() != array.getCapacity())
            return false;
        for (size_t i = 0; i < getLength(); ++i)
            if (operator[](i) != array[i])
                return false;
        return true;
    }

    template<class Derived, class Allocator>
    template<class... Args>
    __host__ __device__ consteval bool ArrayBase<Derived, Allocator>::isTrivialDefaultConstruct() {
        return (sizeof...(Args) == 0) && std::is_trivially_default_constructible<value_type>::value;
    }
}
