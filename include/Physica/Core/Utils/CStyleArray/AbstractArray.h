/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ABSTRACTARRAY_H
#define PHYSICA_ABSTRACTARRAY_H

#include <cstdlib>
#include "QTypeInfo"

namespace Physica::Core {
    enum { Dynamic = 0 };
    /*!
     * \class ArrayArr is the class that handles public part(arr) of specializations of CStyleArray.
     */
    template<class T>
    class AbstractArray {
    protected:
        T* __restrict arr;

        inline AbstractArray();
        inline explicit AbstractArray(size_t capacity);
        inline AbstractArray(AbstractArray<T>&& array) noexcept;

        inline AbstractArray& operator=(AbstractArray<T>&& array) noexcept;

        inline void swap(AbstractArray<T>& array);
    public:
        AbstractArray(const AbstractArray<T>& array) = delete;
        /* Helpers */
        inline void allocate(const T& t, size_t index);
        inline void allocate(T&& t, size_t index);
    };
}

#include "AbstractArrayImpl.h"

#endif