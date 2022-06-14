/*
 * Copyright 2022 WeiBo He.
 *
 * This file is part of Physica.
 *
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

#include "Point.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Lo S. H. Finite Element Mesh Generation[M]. Taylor and Francis, 2015.19-20
     */
    template<size_t dim, class ScalarType>
    class GeoBase {
        using PointType = Point<dim, ScalarType>;
        using VectorType = typename PointType::VectorType;
    public:
        [[nodiscard]] static ScalarType distToSegment(const PointType& p, const PointType& segFrom, const PointType& segTo);
    };

    template<size_t dim, class ScalarType>
    ScalarType GeoBase<dim, ScalarType>::distToSegment(const PointType& p, const PointType& segFrom, const PointType& segTo) {
        VectorType ab = segTo.v() - segFrom.v();
        const ScalarType norm_ab = ab.norm();
        ab *= reciprocal(norm_ab);
        const VectorType ap = p.v() - segFrom.v();
        const ScalarType alpha = ap * ab;

        if (alpha.isNegative())
            return ap.norm();
        if (alpha > norm_ab)
            return (p.v() - segTo.v()).norm();
        return (ap - alpha * ab).norm();
    }
}
