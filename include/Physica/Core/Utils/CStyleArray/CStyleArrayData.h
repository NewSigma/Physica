/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_CSTYLEARRAYDATA_H
#define PHYSICA_CSTYLEARRAYDATA_H

#include <cstdlib>
#include "QTypeInfo"

namespace Physica::Core {
    template<class T>
    class CStyleArrayData {
    protected:
        T* __restrict arr;

        inline CStyleArrayData();
        inline explicit CStyleArrayData(size_t capacity);
        inline CStyleArrayData(CStyleArrayData<T>&& array) noexcept;

        inline CStyleArrayData& operator=(CStyleArrayData<T>&& array) noexcept;

        inline void swap(CStyleArrayData<T>& array);
    public:
        CStyleArrayData(const CStyleArrayData<T>& array) = delete;
        /* Helpers */
        inline void allocate(const T& t, size_t index);
        inline void allocate(T&& t, size_t index);
    };
}

#include "CStyleArrayDataImpl.h"

#endif