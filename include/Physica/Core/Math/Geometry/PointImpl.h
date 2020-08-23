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
#ifndef PHYSICA_POINTIMPL_H
#define PHYSICA_POINTIMPL_H

namespace Physica::Core {
    ////////////////////////////////////////nD////////////////////////////////////////
    template<size_t dim, ScalarType type, bool errorTrack>
    Point<dim, type, errorTrack>::Point(const Point<dim, type, errorTrack>& p) : arr(new Scalar<type, errorTrack>[dim]) {
        for(size_t i = 0; i < length; ++i)
            arr[i] = p.arr[i];
    }

    template<size_t dim, ScalarType type, bool errorTrack>
    Point<dim, type, errorTrack>& Point<dim, type, errorTrack>::operator=(const Point<dim, type, errorTrack>& p) {
        if(&p != this) {
            for(size_t i = 0; i < length; ++i)
                arr[i] = p.arr[i];
        }
        return *this;
    }
    ////////////////////////////////////////1D////////////////////////////////////////
    template<ScalarType type, bool errorTrack>
    Point<1, type, errorTrack>::Point(Scalar<type, errorTrack> x) : x(x) {}

    template<ScalarType type, bool errorTrack>
    Point<1, type, errorTrack>::Point(const Point<1, type, errorTrack>& p) : x(p.x) {}

    template<ScalarType type, bool errorTrack>
    Point<1, type, errorTrack>::Point(Point<1, type, errorTrack>&& p) noexcept : x(std::move(p.x)) {}

    template<ScalarType type, bool errorTrack>
    Point<1, type, errorTrack>& Point<1, type, errorTrack>::operator=(const Point<1, type, errorTrack>& p) {
        x = p.x;
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    Point<1, type, errorTrack>& Point<1, type, errorTrack>::operator=(Point<1, type, errorTrack>&& p) noexcept {
        x = std::move(p.x);
        return *this;
    }
    ////////////////////////////////////////2D////////////////////////////////////////
    template<ScalarType type, bool errorTrack>
    Point<2, type, errorTrack>::Point(Scalar<type, errorTrack> x, Scalar<type, errorTrack> y) : x(x), y(y) {}

    template<ScalarType type, bool errorTrack>
    Point<2, type, errorTrack>::Point(const Point<2, type, errorTrack>& p) : x(p.x), y(p.y) {}

    template<ScalarType type, bool errorTrack>
    Point<2, type, errorTrack>::Point(Point<2, type, errorTrack>&& p) noexcept : x(std::move(p.x)), y(std::move(p.y)) {}

    template<ScalarType type, bool errorTrack>
    Point<2, type, errorTrack>& Point<2, type, errorTrack>::operator=(const Point<2, type, errorTrack>& p) {
        x = p.x;
        y = p.y;
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    Point<2, type, errorTrack>& Point<2, type, errorTrack>::operator=(Point<2, type, errorTrack>&& p) noexcept {
        x = std::move(p.x);
        y = std::move(p.y);
        return *this;
    }
    ////////////////////////////////////////3D////////////////////////////////////////
    template<ScalarType type, bool errorTrack>
    Point<3, type, errorTrack>::Point(Scalar<type, errorTrack> x, Scalar<type, errorTrack> y, Scalar<type, errorTrack> z) : x(x), y(y), z(z) {}

    template<ScalarType type, bool errorTrack>
    Point<3, type, errorTrack>::Point(const Point<3, type, errorTrack>& p) : x(p.x), y(p.y), z(p.z) {}

    template<ScalarType type, bool errorTrack>
    Point<3, type, errorTrack>::Point(Point<3, type, errorTrack>&& p) noexcept : x(std::move(p.x)), y(std::move(p.y)), z(std::move(p.z)) {}

    template<ScalarType type, bool errorTrack>
    Point<3, type, errorTrack>& Point<3, type, errorTrack>::operator=(const Point<3, type, errorTrack>& p) {
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }

    template<ScalarType type, bool errorTrack>
    Point<3, type, errorTrack>& Point<3, type, errorTrack>::operator=(Point<3, type, errorTrack>&& p) noexcept {
        x = std::move(p.x);
        y = std::move(p.y);
        z = std::move(p.z);
        return *this;
    }
}

#endif