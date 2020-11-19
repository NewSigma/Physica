/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_ABSTRACTCSTYLEARRAY_H
#define PHYSICA_ABSTRACTCSTYLEARRAY_H

#include <cstddef>

namespace Physica::Utils {
    /*!
     * Public parts among specializations of CStyleArray.
     */
    template<class T>
    class AbstractCStyleArray {
    public:
        class Iterator {
            T* p;
        public:
            ~Iterator() = default;

            Iterator(const Iterator& ite) : p(ite.p) {}
            Iterator(Iterator&& ite) noexcept : p(ite.p) {}

            /* Operators */
            Iterator& operator=(const Iterator& ite);
            Iterator& operator=(Iterator&& ite) noexcept;
            bool operator==(const Iterator& ite) const noexcept { return p == ite.p; }
            bool operator!=(const Iterator& ite) const noexcept { return p != ite.p; }
            Iterator& operator++();
            const Iterator operator++(int);
            T& operator*() { return *p; }
        private:
            explicit Iterator(T* p) : p(p) {}

            friend class AbstractCStyleArray;
        };
    protected:
        T* __restrict arr;
        size_t length;
    public:
        AbstractCStyleArray() = delete;
        AbstractCStyleArray(const AbstractCStyleArray& array) = delete;
        AbstractCStyleArray(AbstractCStyleArray&& array) noexcept;
        ~AbstractCStyleArray();
        /* Operators */
        AbstractCStyleArray& operator=(const AbstractCStyleArray& array) = delete;
        AbstractCStyleArray& operator=(AbstractCStyleArray&& array) noexcept;
        T& operator[](size_t index) { Q_ASSERT(index < length); return arr[index]; }
        const T& operator[](size_t index) const { Q_ASSERT(index < length); return arr[index]; }
        bool operator==(const AbstractCStyleArray& array) const;
        bool operator!=(const AbstractCStyleArray& array) const { return !(operator==(array)); }
        /* Iterator */
        Iterator begin() { return Iterator(arr); }
        Iterator end() { return Iterator(arr + length); }
        /* Helpers */
        inline void allocate(const T& t, size_t index);
        inline void allocate(T&& t, size_t index);
        T cutLast();
        void swap(AbstractCStyleArray& array);
        /* Getters */
        [[nodiscard]] size_t getLength() const { return length; }
        [[nodiscard]] bool empty() const { return length == 0; }
        /* Setters */
        /*!
         * Low level api. Designed for performance.
         * \size must larger than current length. Because we can not delete the elements we do not need if not.
         * Elements between old length and \size have not allocated. DO NOT try to visit them.
         */
        void setLength(size_t size) { Q_ASSERT(length <= size); length = size; }
    protected:
        /*!
         * @p arr should be allocated by its subclasses using malloc.
         */
        AbstractCStyleArray(T* __restrict arr, size_t length) : arr(arr), length(length) {}
        AbstractCStyleArray(std::initializer_list<T> list);
    };

    template<class T>
    inline void swap(AbstractCStyleArray<T>& a1, AbstractCStyleArray<T>& a2) {
        a1.swap(a2);
    }
}

#include "AbstractCStyleArrayImpl.h"

#endif