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

#include "Physica/Core/Math/Geometry/GeoBase2D.h"

namespace Physica::Core {
    template<class ScalarType>
    Triangle1<ScalarType>::Triangle1(PosArray pos_, IndexArray globalNodes)
            : Base(std::move(globalNodes)), pos(std::move(pos_)) {}

    template<class ScalarType>
    Triangle1<ScalarType>& Triangle1<ScalarType>::operator=(Triangle1<ScalarType> elem) noexcept {
        swap(elem);
        return *this;
    }

    template<class ScalarType>
    typename Triangle1<ScalarType>::MatrixType Triangle1<ScalarType>::jacobi([[maybe_unused]] VectorType localPos) const {
        return MatrixType{pos[1][0] - pos[0][0], pos[1][1] - pos[0][1], pos[2][0] - pos[0][0], pos[2][1] - pos[0][1]};
    }

    template<class ScalarType>
    typename Triangle1<ScalarType>::MatrixType Triangle1<ScalarType>::inv_jacobi([[maybe_unused]] VectorType globalPos) const {
        const VectorType delta21 = pos[1] - pos[0];
        const VectorType delta31 = pos[2] - pos[0];
        const ScalarType factor = reciprocal(delta21[0] * delta31[1] - delta31[0] * delta21[1]);
        return MatrixType{delta31[1], -delta21[1], -delta31[0], delta21[0]} * factor;
    }

    template<class ScalarType>
    bool Triangle1<ScalarType>::contains(const VectorType& point) const {
        return GeoBase2D<ScalarType>::pointOnPoly(point, pos);
    }

    template<class ScalarType>
    typename Triangle1<ScalarType>::VectorType Triangle1<ScalarType>::getNodePos(size_t localNode) const {
        return pos[localNode];
    }

    template<class ScalarType>
    typename Triangle1<ScalarType>::VectorType Triangle1<ScalarType>::toLocalPos(VectorType globalPos) const {
        const VectorType delta21 = pos[1] - pos[0];
        const VectorType delta31 = pos[2] - pos[0];
        globalPos -= pos[0];
        const ScalarType factor = reciprocal(delta21[0] * delta31[1] - delta31[0] * delta21[1]);
        return VectorType{globalPos[0] * delta31[1] - globalPos[1] * delta31[0], globalPos[1] * delta21[0] - globalPos[0] * delta21[1]} * factor;
    }

    template<class ScalarType>
    typename Triangle1<ScalarType>::VectorType Triangle1<ScalarType>::toGlobalPos(VectorType localPos) const {
        return (ScalarType::One() - localPos[0] - localPos[1]) * pos[0] + localPos[0] * pos[1] + localPos[1] * pos[2];
    }

    template<class ScalarType>
    void Triangle1<ScalarType>::swap(Triangle1& elem) noexcept {
        Base::swap(elem);
        pos.swap(elem.pos);
    }

    template<class ScalarType>
    ScalarType Triangle1<ScalarType>::baseFunc(size_t localNode, VectorType p) {
        ScalarType temp;
        switch (localNode) {
            case 0:
                temp = ScalarType::One() - p[0] - p[1];
                break;
            case 1:
                temp = p[0];
                break;
            case 2:
                temp = p[1];
        }
        return temp;
    }

    template<class ScalarType>
    ScalarType Triangle1<ScalarType>::dBase_dr(size_t localNode, [[maybe_unused]] VectorType p){
        int temp;
        switch (localNode) {
            case 0:
                temp = -1;
                break;
            case 1:
                temp = 1;
                break;
            case 2:
                temp = 0;
        }
        return ScalarType(temp);
    }

    template<class ScalarType>
    ScalarType Triangle1<ScalarType>::dBase_ds(size_t localNode, [[maybe_unused]] VectorType p) {
        int temp;
        switch (localNode) {
            case 0:
                temp = -1;
                break;
            case 1:
                temp = 0;
                break;
            case 2:
                temp = 1;
        }
        return ScalarType(temp);
    }

    template<class ScalarType>
    typename Triangle1<ScalarType>::VectorType Triangle1<ScalarType>::grad(size_t localNode, VectorType p) {
        return {dBase_dr(localNode, p), dBase_ds(localNode, p)};
    }

    template<class ScalarType>
    template<class Functor>
    ScalarType Triangle1<ScalarType>::integral(Functor func) {
        constexpr static double factor = 0.577350269189626;
        constexpr static double factor1 = (1 - factor) * 0.25;
        constexpr static double factor2 = (1 + factor) * 0.5;
        constexpr static double factor3 = (1 + factor) * 0.25;
        constexpr static double factor4 = (1 - factor) * 0.5;
        return func(VectorType{factor1, factor2}) * factor1 + func(VectorType{factor3, factor4}) * factor3;
    }

    template<class ScalarType>
    Mesh<Triangle1<ScalarType>> Triangle1<ScalarType>::rectangle(VectorType bottomLeft,
                                                                 VectorType topRight,
                                                                 size_t numSeparateX,
                                                                 size_t numSeparateY) {
        const size_t numNodeX = numSeparateX + 1;
        const size_t numNodeY = numSeparateY + 1;
        Mesh<Triangle1<ScalarType>> mesh(numSeparateX * numSeparateY * 2, numNodeX * numNodeY);
        const ScalarType xPerElem = (topRight[0] - bottomLeft[0]) / ScalarType(numSeparateX);
        const ScalarType yPerElem = (topRight[1] - bottomLeft[1]) / ScalarType(numSeparateY);

        VectorType p = bottomLeft;
        size_t nextElem = 0;
        for (size_t y = 0; y < numSeparateY; ++y) {
            for (size_t x = 0; x < numSeparateX; ++x) {
                const size_t nodeBottomLeft = numNodeX * y + x;
                const size_t nodeBottomRight = nodeBottomLeft + 1;
                const size_t nodeTopLeft = nodeBottomLeft + numNodeX;
                const size_t nodeTopRight = nodeBottomRight + numNodeX;
                const VectorType posBottomRight{p[0] + xPerElem, p[1]};
                const VectorType posTopLeft{p[0], p[1] + yPerElem};
                const VectorType posTopRight{p[0] + xPerElem, p[1] + yPerElem};
                mesh.setElem(Triangle1({p, posBottomRight, posTopLeft}
                                      , {nodeBottomLeft, nodeBottomRight, nodeTopLeft})
                            , nextElem++);
                mesh.setElem(Triangle1({posBottomRight, posTopRight, posTopLeft}
                                      , {nodeBottomRight, nodeTopRight, nodeTopLeft})
                            , nextElem++);
                p[0] += xPerElem;
            }
            p[1] += yPerElem;
            p[0] = bottomLeft[0];
        }
        return mesh;
    }
}
