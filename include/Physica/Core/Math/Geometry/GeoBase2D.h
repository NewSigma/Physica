/*
 * Copyright 2022-2024 Weibo He.
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

#include "GeoBase.h"

namespace Physica {
    template<Scalar T>
    class GeoBase2D : public GeoBase<2, T> {
    public:
        using Base = GeoBase<2, T>;
        using typename Base::PointType;
        using typename Base::VectorType;
    public:
        template<size_t PolyVertex>
        [[nodiscard]] static bool pointInPoly(const VectorType& p, const Array<VectorType, PolyVertex>& vertexes);
        template<size_t PolyVertex>
        [[nodiscard]] static bool pointOnPoly(const VectorType& p, const Array<VectorType, PolyVertex>& vertexes);
    private:
        template<size_t PolyVertex>
        static T pointPolyImpl(const VectorType& p, const Array<VectorType, PolyVertex>& vertexes);
    };
    /**
     * \param vertexes
     * The vertexes are sorted anticlock anticlockwise.
     *
     * Reference:
     * [1] Lo S. H. Finite Element Mesh Generation[M]. Taylor and Francis, 2015.32
     */
    template<Scalar T>
    template<size_t PolyVertex>
    bool GeoBase2D<T>::pointInPoly(const VectorType& p, const Array<VectorType, PolyVertex>& vertexes) {
        const size_t numVertex = vertexes.getLength();
        for (size_t i = 0; i < numVertex; ++i) {
            const size_t next_i = (i + 1) % numVertex;
            const VectorType ab = vertexes[next_i] - vertexes[i];
            const VectorType bp = p - vertexes[next_i];
            if (ab[0] * bp[1] <= ab[1] * bp[0])
                return false;
        }
        return true;
    }

    template<Scalar T>
    template<size_t PolyVertex>
    bool GeoBase2D<T>::pointOnPoly(const VectorType& p, const Array<VectorType, PolyVertex>& vertexes) {
        const size_t numVertex = vertexes.getLength();
        for (size_t i = 0; i < numVertex; ++i) {
            const size_t next_i = (i + 1) % numVertex;
            const VectorType ab = vertexes[next_i] - vertexes[i];
            const VectorType bp = p - vertexes[next_i];
            if (ab[0] * bp[1] < ab[1] * bp[0])
                return false;
        }
        return true;
    }
}
