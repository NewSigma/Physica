/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POINT_H
#define PHYSICA_POINT_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    template<size_t dim = 3, ScalarType type = MultiPrecision, bool errorTrack = true>
    class Point {
        Scalar<type, errorTrack>* arr;
    public:
        static constexpr size_t length = dim;

        Point() : arr(new Scalar<type, errorTrack>[dim]) {}
        Point(const Point& p);
        Point(Point&& p) noexcept : arr(p.arr) { p.arr = nullptr; }
        ~Point() { delete[] arr; }
        /* Operators */
        Point& operator=(const Point& p);
        Point& operator=(Point&& p) noexcept { this->~Point(); arr = p.arr; p.arr = nullptr; }
    };

    template<ScalarType type, bool errorTrack>
    class Point<1, type, errorTrack> {
    public:
        Scalar<type, errorTrack> x;
    };

    template<ScalarType type, bool errorTrack>
    class Point<2, type, errorTrack> {
        Scalar<type, errorTrack> x, y;
    };

    template<ScalarType type, bool errorTrack>
    class Point<3, type, errorTrack> {
        Scalar<type, errorTrack> x, y, z;
    };
}

#include "PointImpl.h"

#endif