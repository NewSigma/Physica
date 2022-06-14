/*
 * Copyright 2020-2022 WeiBo He.
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
        using VectorType = Vector<ScalarType, dim>;
    private:
        VectorType vec;
    public:
        Point() = default;
        template<class Derived>
        Point(const LValueVector<Derived>& vec_) : vec(vec_) {}
        Point(std::initializer_list<ScalarType> list) : vec(std::move(list)) {}
        Point(const Point& p) = default;
        Point(Point&& p) noexcept = default;
        ~Point() = default;
        /* Operators */
        Point& operator=(const Point& p) = default;
        Point& operator=(Point&& p) noexcept = default;
        /* Getters */
        VectorType& v() noexcept { return vec; }
        const VectorType& v() const noexcept { return vec; }
        const ScalarType& x() const { return vec[0]; }
        const ScalarType& y() const { return vec[1]; }
        const ScalarType& z() const { return vec[2]; }
        const ScalarType& w() const { return vec[3]; }
        ScalarType dist(const Point& p) const;
    };

    template<size_t dim, class ScalarType>
    std::ostream& operator<<(std::ostream& os, const Point<dim, ScalarType>& p) {
        return os << p.v().format().setPrefix("(").setSeparator(", ").setSuffix(")");
    }

    template<size_t dim, class ScalarType>
    ScalarType Point<dim, ScalarType>::dist(const Point<dim, ScalarType>& p) const {
        return (vec - p.vec).norm();
    }

    typedef Point<1> Point1D;
    typedef Point<2> Point2D;
    typedef Point<3> Point3D;
}
