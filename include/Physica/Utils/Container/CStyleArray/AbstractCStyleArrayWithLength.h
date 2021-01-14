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

#include "AbstractCStyleArray.h"

namespace Physica::Utils::Intenal {
    template<class Derived>
    class AbstractCStyleArrayWithLength : public AbstractCStyleArray<Derived> {
    private:
        using Base = AbstractCStyleArray<Derived>;
        using typename Base::T;
    protected:
        using Base::arr;
        size_t length;
    public:
        AbstractCStyleArrayWithLength() = delete;
        ~AbstractCStyleArrayWithLength() = default;
        /* Operators */
        Derived& operator<<(const T& t) { Base::getDerived().append(t); return Base::getDerived(); }
        Derived& operator<<(T&& t) { Base::getDerived().append(std::move(t)); return Base::getDerived(); }
        Derived& operator<<(const Derived& array) { Base::getDerived().append(array); return Base::getDerived(); }
        Derived& operator<<(Derived&& array) { Base::getDerived().append(std::move(array)); return Base::getDerived(); }
        /* Helpers */
        T cutLast();
        inline void grow(const T& t);
        inline void grow(T&& t);
        void removeAt(size_t index);
        /* Setters */
        /**
         * Low level api. Designed for performance.
         * \size must larger than current length. Because we can not delete the elements we do not need if not.
         * Elements between old length and \size have not allocated. DO NOT try to visit them.
         */
        void setLength(size_t size) { Q_ASSERT(length <= size); length = size; }
    protected:
        AbstractCStyleArrayWithLength(size_t length_, size_t capacity);
        AbstractCStyleArrayWithLength(const AbstractCStyleArrayWithLength& array);
        AbstractCStyleArrayWithLength(AbstractCStyleArrayWithLength&& array) noexcept;
        /* Operators */
        AbstractCStyleArrayWithLength& operator=(const AbstractCStyleArrayWithLength& array);
        AbstractCStyleArrayWithLength& operator=(AbstractCStyleArrayWithLength&& array) noexcept;
        /* Helpers */
        void swap(AbstractCStyleArrayWithLength& array);
    };
}

#include "AbstractCStyleArrayWithLengthImpl.h"