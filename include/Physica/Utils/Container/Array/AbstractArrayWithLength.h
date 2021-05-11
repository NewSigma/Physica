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

#include "AbstractArray.h"

namespace Physica::Utils::Internal {
    template<class Derived>
    class AbstractArrayWithLength : public AbstractArray<Derived> {
    private:
        using Base = AbstractArray<Derived>;
        using typename Base::T;
    protected:
        using Base::arr;
        size_t length;
    public:
        AbstractArrayWithLength() = delete;
        ~AbstractArrayWithLength() = default;
        /* Operators */
        AbstractArrayWithLength& operator=(const AbstractArrayWithLength& array) = delete;
        AbstractArrayWithLength& operator=(AbstractArrayWithLength&& array) noexcept = delete;
        Derived& operator<<(const T& t) { Base::getDerived().append(t); return Base::getDerived(); }
        Derived& operator<<(T&& t) { Base::getDerived().append(std::move(t)); return Base::getDerived(); }
        Derived& operator<<(const Derived& array) { Base::getDerived().append(array); return Base::getDerived(); }
        Derived& operator<<(Derived&& array) { Base::getDerived().append(std::move(array)); return Base::getDerived(); }
        /* Operations */
        using Base::allocate;
        T cutLast();
        inline void grow(const T& t);
        inline void grow(T&& t);
        void removeAt(size_t index);
        void clear() noexcept;
        void insert(const T& t, size_t index);
        void insert(T&& t, size_t index);
        /* Setters */
        /**
         * Low level api. Designed for performance.
         * \size must larger than current length. Because we can not delete the elements we do not need if not.
         * Elements between old length and \size have not allocated. DO NOT try to visit them.
         */
        void setLength(size_t size) { assert(length <= size && size <= Base::getDerived().getCapacity()); length = size; }
    protected:
        explicit AbstractArrayWithLength(size_t capacity);
        AbstractArrayWithLength(size_t length_, size_t capacity);
        AbstractArrayWithLength(size_t length_, T* __restrict arr_);
        AbstractArrayWithLength(const AbstractArrayWithLength& array);
        AbstractArrayWithLength(AbstractArrayWithLength&& array) noexcept;
        /* Helpers */
        void swap(AbstractArrayWithLength& array);
    };
}

#include "AbstractArrayWithLengthImpl.h"