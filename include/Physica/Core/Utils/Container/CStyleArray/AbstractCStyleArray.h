/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ABSTRACTCSTYLEARRAY_H
#define PHYSICA_ABSTRACTCSTYLEARRAY_H

#include <cstddef>

namespace Physica::Core {
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
