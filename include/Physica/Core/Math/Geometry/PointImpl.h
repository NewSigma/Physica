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
    template<size_t dim, ScalarOption option, bool errorTrack>
    Point<dim, option, errorTrack>::Point(const Point<dim, option, errorTrack>& p) : arr(new Scalar<option, errorTrack>[dim]) {
        for(size_t i = 0; i < length; ++i)
            arr[i] = p.arr[i];
    }

    template<size_t dim, ScalarOption option, bool errorTrack>
    Point<dim, option, errorTrack>& Point<dim, option, errorTrack>::operator=(const Point<dim, option, errorTrack>& p) {
        if(&p != this) {
            for(size_t i = 0; i < length; ++i)
                arr[i] = p.arr[i];
        }
        return *this;
    }
    ////////////////////////////////////////1D////////////////////////////////////////
    template<ScalarOption option, bool errorTrack>
    Point<1, option, errorTrack>::Point(Scalar<option, errorTrack> x) : x(x) {}

    template<ScalarOption option, bool errorTrack>
    Point<1, option, errorTrack>::Point(const Point<1, option, errorTrack>& p) : x(p.x) {}

    template<ScalarOption option, bool errorTrack>
    Point<1, option, errorTrack>::Point(Point<1, option, errorTrack>&& p) noexcept : x(std::move(p.x)) {}

    template<ScalarOption option, bool errorTrack>
    Point<1, option, errorTrack>& Point<1, option, errorTrack>::operator=(const Point<1, option, errorTrack>& p) {
        x = p.x;
        return *this;
    }

    template<ScalarOption option, bool errorTrack>
    Point<1, option, errorTrack>& Point<1, option, errorTrack>::operator=(Point<1, option, errorTrack>&& p) noexcept {
        x = std::move(p.x);
        return *this;
    }
    ////////////////////////////////////////2D////////////////////////////////////////
    template<ScalarOption option, bool errorTrack>
    Point<2, option, errorTrack>::Point(Scalar<option, errorTrack> x, Scalar<option, errorTrack> y) : x(x), y(y) {}

    template<ScalarOption option, bool errorTrack>
    Point<2, option, errorTrack>::Point(const Point<2, option, errorTrack>& p) : x(p.x), y(p.y) {}

    template<ScalarOption option, bool errorTrack>
    Point<2, option, errorTrack>::Point(Point<2, option, errorTrack>&& p) noexcept : x(std::move(p.x)), y(std::move(p.y)) {}

    template<ScalarOption option, bool errorTrack>
    Point<2, option, errorTrack>& Point<2, option, errorTrack>::operator=(const Point<2, option, errorTrack>& p) {
        x = p.x;
        y = p.y;
        return *this;
    }

    template<ScalarOption option, bool errorTrack>
    Point<2, option, errorTrack>& Point<2, option, errorTrack>::operator=(Point<2, option, errorTrack>&& p) noexcept {
        x = std::move(p.x);
        y = std::move(p.y);
        return *this;
    }
    ////////////////////////////////////////3D////////////////////////////////////////
    template<ScalarOption option, bool errorTrack>
    Point<3, option, errorTrack>::Point(Scalar<option, errorTrack> x, Scalar<option, errorTrack> y, Scalar<option, errorTrack> z) : x(x), y(y), z(z) {}

    template<ScalarOption option, bool errorTrack>
    Point<3, option, errorTrack>::Point(const Point<3, option, errorTrack>& p) : x(p.x), y(p.y), z(p.z) {}

    template<ScalarOption option, bool errorTrack>
    Point<3, option, errorTrack>::Point(Point<3, option, errorTrack>&& p) noexcept : x(std::move(p.x)), y(std::move(p.y)), z(std::move(p.z)) {}

    template<ScalarOption option, bool errorTrack>
    Point<3, option, errorTrack>& Point<3, option, errorTrack>::operator=(const Point<3, option, errorTrack>& p) {
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }

    template<ScalarOption option, bool errorTrack>
    Point<3, option, errorTrack>& Point<3, option, errorTrack>::operator=(Point<3, option, errorTrack>&& p) noexcept {
        x = std::move(p.x);
        y = std::move(p.y);
        z = std::move(p.z);
        return *this;
    }
}

#endif