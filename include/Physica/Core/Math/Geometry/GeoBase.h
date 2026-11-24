/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Point.h"

namespace Physica {
    /**
     * Reference:
     * [1] Lo S. H. Finite Element Mesh Generation[M]. Taylor and Francis, 2015.19-20
     */
    template<size_t Dim, Scalar T>
    class GeoBase {
    public:
        using PointType = Point<Dim, T>;
        using VectorType = PointType::VectorType;
    public:
        [[nodiscard]] static T distToSegment(const PointType& p, const PointType& segFrom, const PointType& segTo);
    };

    template<size_t Dim, Scalar T>
    T GeoBase<Dim, T>::distToSegment(const PointType& p, const PointType& segFrom, const PointType& segTo) {
        VectorType ab = segTo.v() - segFrom.v();
        const T norm_ab = ab.norm();
        ab *= reciprocal(norm_ab);
        const VectorType ap = p.v() - segFrom.v();
        const T alpha = ap * ab;

        if (alpha.isNegative())
            return ap.norm();
        if (alpha > norm_ab)
            return (p.v() - segTo.v()).norm();
        return (ap - alpha * ab).norm();
    }
}
