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

#include <cstring>
#include "ArrayBase.h"

namespace Physica {
    template<class Derived, class Allocator>
    __host__ __device__ inline auto ArrayBase<Derived, Allocator>::operator[](size_t index) -> lvalue_reference {
        assert(index < getLength() && "[Error]: Index overflow");
        return data()[index];
    }

    template<class Derived, class Allocator>
    __host__ __device__ inline auto ArrayBase<Derived, Allocator>::operator[](size_t index) const -> const_lvalue_reference {
        assert(index < getLength() && "[Error]: Index overflow");
        return data()[index];
    }

    template<class Derived, class Allocator>
    __host__ __device__ bool ArrayBase<Derived, Allocator>::operator==(const ArrayBase& array) const {
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
    void ArrayBase<Derived, Allocator>::read(const auto& loc, const char* name) {
        const auto group = loc.openGroup(name);
        char buffer[32]; //32 is enough for uint64_t
        for (size_t i = 0; i < getLength(); ++i) {
            std::sprintf(buffer, "%zu", i);
            (*this)[i].read(group, buffer);
        }
    }

    template<class Derived, class Allocator>
    void ArrayBase<Derived, Allocator>::write(auto& loc, const char* name) const {
        auto group = loc.openGroup(name);
        char buffer[32]; //32 is enough for uint64_t
        for (size_t i = 0; i < getLength(); ++i) {
            std::sprintf(buffer, "%zu", i);
            (*this)[i].write(group, buffer);
        }
    }

    template<class Derived, class Allocator>
    __host__ __device__ auto ArrayBase<Derived, Allocator>::data_ptr(size_t index) noexcept -> pointer {
        assert(index < getLength());
        return data() + index;
    }

    template<class Derived, class Allocator>
    __host__ __device__ auto ArrayBase<Derived, Allocator>::data_ptr(size_t index) const noexcept -> const_pointer {
        return const_cast<This&>(*this).data_ptr(index);
    }

    template<class Derived, class Allocator>
    template<class... Args>
    __host__ __device__ consteval bool ArrayBase<Derived, Allocator>::isTrivialDefaultConstruct() {
        return (sizeof...(Args) == 0) && std::is_trivially_default_constructible<value_type>::value;
    }
}
