/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POINTIMPL_H
#define PHYSICA_POINTIMPL_H

namespace Physica::Core {
    template<size_t dim, ScalarType type, bool errorTrack>
    Point<dim, type, errorTrack>(const Point<dim, type, errorTrack>& p) : arr(new Scalar<type, errorTrack>[dim]) {
        for(size_t i = 0; i < length; ++i)
            arr[i] = p.arr[i];
    }

    template<size_t dim, ScalarType type, bool errorTrack>
    Point<dim, type, errorTrack>& operator=(const Point<dim, type, errorTrack>& p) {
        if(&p != this) {
            for(size_t i = 0; i < length; ++i)
                arr[i] = p.arr[i];
        }
        return *this;
    }
}

#endif