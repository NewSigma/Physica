/*
 * Copyright 2023 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    /**
     * Calculate cross section area of plain ax + by + cz + d = 0 and cube [-1, 1] * [-1, 1] * [-1, 1]
     * Algorithm introduced in [1] and its erratum [2]
     * 
     * Reference:
     * [1] Phys. Rev. 144, 390 (1966); https://doi.org/10.1103/PhysRev.144.390
     * [2] Phys. Rev. 147, 670 (1966); https://doi.org/10.1103/PhysRev.147.670.2
     */
    template<class ScalarType>
    class CubeCross {
        using Vector3D = Vector<ScalarType, 3>;
        using Vector4D = Vector<ScalarType, 4>;
    public:
        enum CrossType : char {
            Parallelogram,
            Hexagon,
            Pentagon,
            Quadrangle,
            Triangle,
            None
        };
    private:
        ScalarType area;
        CrossType type;
    public:
        CubeCross(Vector4D coeff);
        CubeCross(const CubeCross&) = default;
        CubeCross(CubeCross&&) noexcept = default;
        ~CubeCross() = default;
        /* Operators */
        CubeCross& operator=(CubeCross obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void swap(CubeCross& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] ScalarType getArea() const noexcept { return area; }
        [[nodiscard]] CrossType getType() const noexcept { return type; }
    private:
        static Vector3D sortVector(Vector3D v);
    };

    template<class ScalarType>
    CubeCross<ScalarType>::CubeCross(Vector4D coeff) {
        assert(!coeff.squaredNorm().isZero() && "[Error]: Invalid plain");
        const Vector3D coeff1 = coeff.head(3);
        const ScalarType repNorm = reciprocal(coeff1.norm());
        const ScalarType dist = abs(coeff[3]) * repNorm;
        const Vector3D normalVec = sortVector(abs(coeff1) * repNorm);

        const ScalarType signedDist1 = normalVec[0] - normalVec[1] - normalVec[2];
        const ScalarType dist1 = abs(signedDist1);
        const ScalarType dist2 = normalVec[0] - normalVec[1] + normalVec[2];
        const ScalarType dist3 = normalVec[0] + normalVec[1] - normalVec[2];
        const ScalarType dist4 = normalVec.sum();
        if (dist > dist4) {
            area = 0;
            type = None;
            return;
        }

        const ScalarType factor = normalVec[0] * normalVec[1] * normalVec[2];
        if (dist < dist2) {
            if (dist < dist1) {
                if (!signedDist1.isNegative()) {
                    area = ScalarType(4) / normalVec[0];
                    type = Parallelogram;
                }
                else {
                    area = normalVec[0] * normalVec[1] + normalVec[0] * normalVec[2] + normalVec[1] * normalVec[2];
                    area = ScalarType(2) * area - (square(dist) + ScalarType(1));
                    area /= factor;
                    type = Hexagon;
                }
            }
            else {
                area = normalVec[0] * normalVec[1] + normalVec[0] * normalVec[2] + normalVec[1] * normalVec[2] * ScalarType(3);
                area += dist * signedDist1 - (square(dist) + ScalarType(1)) * ScalarType(0.5);
                area /= factor;
                type = Pentagon;
            }
        }
        else {
            if (dist < dist3) {
                area = normalVec[0] + normalVec[1] - dist;
                area *= ScalarType(2) / (normalVec[0] * normalVec[1]);
                type = Quadrangle;
            }
            else {
                area = square(dist4 - dist) / (factor * ScalarType(2));
                type = Triangle;
            }
        }
    }

    template<class ScalarType>
    void CubeCross<ScalarType>::swap(CubeCross& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        area.swap(obj.area);
        std::swap(type, obj.type);
    }

    template<class ScalarType>
    typename CubeCross<ScalarType>::Vector3D CubeCross<ScalarType>::sortVector(Vector3D v) {
        if (v[0] < v[1])
            v[0].swap(v[1]);
        if (v[0] < v[2])
            v[0].swap(v[2]);
        if (v[1] < v[2])
            v[1].swap(v[2]);
        return v;
    }
}
