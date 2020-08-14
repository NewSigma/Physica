/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_POINT_H
#define PHYSICA_POINT_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    /*!
     * By default, point do not need high precision, so @param type is set to @enum ScalarType::Float.
     */
    template<size_t dim, ScalarType type = Float, bool errorTrack = false>
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
        Scalar<type, errorTrack> x;
    public:
        explicit Point(Scalar<type, errorTrack> x);
        Point(const Point& p);
        Point(Point&& p) noexcept;
        ~Point() = default;
        /* Operators */
        Point& operator=(const Point& p);
        Point& operator=(Point&& p) noexcept;
    };

    template<ScalarType type, bool errorTrack>
    class Point<2, type, errorTrack> {
        Scalar<type, errorTrack> x, y;
    public:
        Point(Scalar<type, errorTrack> x, Scalar<type, errorTrack> y);
        Point(const Point& p);
        Point(Point&& p) noexcept;
        ~Point() = default;
        /* Operators */
        Point& operator=(const Point& p);
        Point& operator=(Point&& p) noexcept;
    };

    template<ScalarType type, bool errorTrack>
    class Point<3, type, errorTrack> {
        Scalar<type, errorTrack> x, y, z;
    public:
        Point(Scalar<type, errorTrack> x, Scalar<type, errorTrack> y, Scalar<type, errorTrack> z);
        Point(const Point& p);
        Point(Point&& p) noexcept;
        ~Point() = default;
        /* Operators */
        Point& operator=(const Point& p);
        Point& operator=(Point&& p) noexcept;
    };

    typedef Point<1> Point1D;
    typedef Point<2> Point2D;
    typedef Point<3> Point3D;
}

#include "PointImpl.h"

#endif