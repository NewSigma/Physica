/*
 * Copyright 2021-2024 Weibo He.
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

#include "ArrayBase.h"

namespace Physica::Utils {
    template<class Derived, class Allocator>
    class DynamicArrayBase : public ArrayBase<Derived, Allocator> {
    public:
        using Base = ArrayBase<Derived, Allocator>;
        using typename Base::allocator_type;
        using typename Base::ValueType;
        using typename Base::pointer;
        using typename Base::const_pointer;
        using typename Base::lvalue_reference;
        using typename Base::const_lvalue_reference;
        using typename Base::rvalue_reference;
    protected:
        pointer __restrict arr;
        allocator_type alloc;
        size_t length;
    public:
        DynamicArrayBase() = delete;
        ~DynamicArrayBase();
        /* Operators */
        Derived& operator<<(const_lvalue_reference t) { Base::getDerived().append(t); return Base::getDerived(); }
        Derived& operator<<(rvalue_reference t) { Base::getDerived().append(std::move(t)); return Base::getDerived(); }
        Derived& operator<<(const Derived& array) { Base::getDerived().append(array); return Base::getDerived(); }
        Derived& operator<<(Derived&& array) { Base::getDerived().append(std::move(array)); return Base::getDerived(); }
        /* Operations */
        ValueType cutLast();
        template<class... Args>
        inline void grow(Args&&... args);
        void removeAt(size_t index);
        void clear() noexcept;
        /* Setters */
        inline void setLength(size_t size);
        /* Getters */
        [[nodiscard]] __host__ __device__ pointer data() noexcept { return arr; }
        [[nodiscard]] __host__ __device__ const_pointer data() const noexcept { return arr; }
        [[nodiscard]] allocator_type get_allocator() const noexcept { return alloc; }
    protected:
        explicit DynamicArrayBase(size_t capacity);
        DynamicArrayBase(size_t length_, size_t capacity);
        DynamicArrayBase(size_t length_, pointer arr_);
        DynamicArrayBase(const DynamicArrayBase& array);
        DynamicArrayBase(DynamicArrayBase&& array) noexcept;
        /* Operators */
        DynamicArrayBase& operator=(DynamicArrayBase array) noexcept;
        /* Operations */
        void swap(DynamicArrayBase& __restrict array) noexcept;
        /* Setters */
        [[nodiscard]] inline pointer release() noexcept;
    };
}

#include "DynamicArrayBaseImpl.h"