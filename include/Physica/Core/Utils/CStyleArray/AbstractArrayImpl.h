/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ABSTRACTARRAYIMPL_H
#define PHYSICA_ABSTRACTARRAYIMPL_H

#ifndef PHYSICA_ABSTRACTARRAY_H
    #include "Core/Header/Utils/CStyleArray/CStyleArrayArr.h"
#endif

namespace Physica::Core {
    template<class T>
    inline AbstractArray<T>::AbstractArray() : arr(nullptr) {}

    template<class T>
    inline AbstractArray<T>::AbstractArray(size_t capacity) : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))) {}

    template<class T>
    inline AbstractArray<T>::AbstractArray(AbstractArray<T> &&array) noexcept : arr(array.arr) {
        array.arr = nullptr;
    }

    template<class T>
    inline AbstractArray<T>& AbstractArray<T>::operator=(AbstractArray<T>&& array) noexcept {
        arr = array.arr;
        array.arr = nullptr;
    }
    /*!
     * Low level api. Designed for performance.
     * Simply allocate a \T at position \index.
     * You must ensure position \index is not used or a memory leak will occur.
     */
    template<class T>
    inline void AbstractArray<T>::allocate(const T& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(t);
        else
            *(arr + index) = t;
    }

    template<class T>
    inline void AbstractArray<T>::allocate(T&& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(std::move(t));
        else
            *(arr + index) = t;
    }

    template<class T>
    inline void AbstractArray<T>::swap(AbstractArray<T>& array) {
        std::swap(arr, array.arr);
    }
}

#endif