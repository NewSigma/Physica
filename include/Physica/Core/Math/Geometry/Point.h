/*
 * Copyright 2020-2021 WeiBo He.
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
#pragma once

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    /*!
     * By default, point do not need high precision, so @param type is set to @enum ScalarType::Float.
     */
    template<size_t dim, class ScalarType = Scalar<Float, false>>
    class Point {
        static_assert(dim > 0, "0 dim point is not allowed\n");
    public:
        static constexpr size_t length = dim;
    private:
        Vector<ScalarType, dim> v;
    public:
        Point() = default;
        template<class VectorType>
        Point(const LValueVector<VectorType>& v_) : v(v_) {}
        Point(std::initializer_list<ScalarType> list) : v(std::move(list)) {}
        Point(const Point& p) = default;
        Point(Point&& p) noexcept = default;
        ~Point() = default;
        /* Operators */
        Point& operator=(const Point& p) = default;
        Point& operator=(Point&& p) noexcept = default;
        /* Getters */
        const ScalarType& getX() const { return v[0]; }
        const ScalarType& getY() const { return v[1]; }
        const ScalarType& getZ() const { return v[2]; }
        const ScalarType& getW() const { return v[3]; }
    };

    typedef Point<1> Point1D;
    typedef Point<2> Point2D;
    typedef Point<3> Point3D;
}
