/*
 * Copyright 2021 WeiBo He.
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

#include "ArrayStorage.h"

namespace Physica::Utils::Internal {
    template<class Derived>
    class DynamicArrayBase : public ArrayStorage<Derived> {
    public:
        using Base = ArrayStorage<Derived>;
        using typename Base::ValueType;
        using typename Base::PointerType;
        using typename Base::LValueReferenceType;
        using typename Base::ConstLValueReferenceType;
        using typename Base::RValueReferenceType;
    protected:
        using Base::arr;
        using Base::alloc;
        size_t length;
    public:
        DynamicArrayBase() = delete;
        ~DynamicArrayBase() = default;
        /* Operators */
        DynamicArrayBase& operator=(const DynamicArrayBase& array) = delete;
        DynamicArrayBase& operator=(DynamicArrayBase&& array) noexcept = delete;
        Derived& operator<<(ConstLValueReferenceType t) { Base::getDerived().append(t); return Base::getDerived(); }
        Derived& operator<<(RValueReferenceType t) { Base::getDerived().append(std::move(t)); return Base::getDerived(); }
        Derived& operator<<(const Derived& array) { Base::getDerived().append(array); return Base::getDerived(); }
        Derived& operator<<(Derived&& array) { Base::getDerived().append(std::move(array)); return Base::getDerived(); }
        /* Operations */
        ValueType cutLast();
        __host__ __device__ inline void grow(ConstLValueReferenceType t);
        __host__ __device__ inline void grow(RValueReferenceType t);
        void removeAt(size_t index);
        __host__ __device__ void clear() noexcept;
        void insert(ConstLValueReferenceType t, size_t index);
        void insert(RValueReferenceType t, size_t index);
        /* Setters */
        /**
         * Low level api. Designed for performance.
         * \size must larger than current length. Because we can not delete the elements we do not need if not.
         * Elements between old length and \size have not allocated. DO NOT try to visit them.
         */
        __host__ __device__ void setLength(size_t size) { assert(length <= size && size <= Base::getDerived().getCapacity()); length = size; }
    protected:
        __host__ __device__ explicit DynamicArrayBase(size_t capacity);
        __host__ __device__ DynamicArrayBase(size_t length_, size_t capacity);
        __host__ __device__ DynamicArrayBase(size_t length_, PointerType arr_);
        __host__ __device__ DynamicArrayBase(const DynamicArrayBase& array);
        __host__ __device__ DynamicArrayBase(DynamicArrayBase&& array) noexcept;
        /* Helpers */
        __host__ __device__ void swap(DynamicArrayBase& array);
    };
}

#include "DynamicArrayBaseImpl.h"