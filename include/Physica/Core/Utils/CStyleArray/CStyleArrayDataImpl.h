/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CSTYLEARRAYDATAIMPL_H
#define PHYSICA_CSTYLEARRAYDATAIMPL_H

#ifndef PHYSICA_CSTYLEARRAYDATA_H
    #include "Core/Header/Utils/CStyleArray/CStyleArrayArr.h"
#endif

namespace Physica::Core {
    template<class T>
    inline CStyleArrayData<T>::CStyleArrayData() : arr(nullptr) {}

    template<class T>
    inline CStyleArrayData<T>::CStyleArrayData(size_t capacity)
            : arr(reinterpret_cast<T*>(malloc(capacity * sizeof(T)))) {}

    template<class T>
    inline CStyleArrayData<T>::CStyleArrayData(CStyleArrayData<T> &&array) noexcept : arr(array.arr) {
        array.arr = nullptr;
    }

    template<class T>
    inline CStyleArrayData<T>& CStyleArrayData<T>::operator=(CStyleArrayData<T>&& array) noexcept {
        arr = array.arr;
        array.arr = nullptr;
    }
    /*!
     * Low level api. Designed for performance.
     * Simply allocate a \T at position \index.
     * You must ensure position \index is not used or a memory leak will occur.
     */
    template<class T>
    inline void CStyleArrayData<T>::allocate(const T& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(t);
        else
            *(arr + index) = t;
    }

    template<class T>
    inline void CStyleArrayData<T>::allocate(T&& t, size_t index) {
        if(QTypeInfo<T>::isComplex)
            new (arr + index) T(std::move(t));
        else
            *(arr + index) = t;
    }

    template<class T>
    inline void CStyleArrayData<T>::swap(CStyleArrayData<T>& array) {
        std::swap(arr, array.arr);
    }
}

#endif