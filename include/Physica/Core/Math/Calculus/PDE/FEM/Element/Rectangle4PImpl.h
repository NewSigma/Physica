/*
 * Copyright 2021-2022 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType>
    Rectangle4P<ScalarType>::Rectangle4P(VectorType bottomLeft_,
                                         VectorType topRight_,
                                         size_t nodeBottomLeft,
                                         size_t nodeBottomRight,
                                         size_t nodeTopLeft,
                                         size_t nodeTopRight)
            : bottomLeft(bottomLeft_), topRight(topRight_) {
        Base::nodes[0] = nodeBottomLeft;
        Base::nodes[1] = nodeBottomRight;
        Base::nodes[2] = nodeTopLeft;
        Base::nodes[3] = nodeTopRight;
    }

    template<class ScalarType>
    Rectangle4P<ScalarType>& Rectangle4P<ScalarType>::operator=(Rectangle4P<ScalarType> elem) noexcept {
        swap(elem);
        return *this;
    }

    template<class ScalarType>
    typename Rectangle4P<ScalarType>::MatrixType Rectangle4P<ScalarType>::jacobi([[maybe_unused]] VectorType localPos) const {
        return MatrixType{(topRight[0] - bottomLeft[0]) * 0.5, 0, 0, (topRight[1] - bottomLeft[1]) * 0.5};
    }

    template<class ScalarType>
    typename Rectangle4P<ScalarType>::MatrixType Rectangle4P<ScalarType>::inv_jacobi([[maybe_unused]] VectorType globalPos) const {
        return MatrixType{ScalarType::Two() / (topRight[0] - bottomLeft[0]), 0, 0, ScalarType::Two() / (topRight[1] - bottomLeft[1])};
    }

    template<class ScalarType>
    bool Rectangle4P<ScalarType>::contains(const VectorType& point) const {
        return bottomLeft[0] <= point[0]
            && point[0] <= topRight[0]
            && bottomLeft[1] <= point[1]
            && point[1] <= topRight[1];
    }

    template<class ScalarType>
    typename Rectangle4P<ScalarType>::VectorType Rectangle4P<ScalarType>::getNodePos(size_t localNode) const {
        switch (localNode) {
            case BottomLeft:
                return bottomLeft;
            case BottomRight:
                return {topRight[0], bottomLeft[1]};
            case TopLeft:
                return {bottomLeft[0], topRight[1]};
            case TopRight:
                return topRight;
        }
        throw std::invalid_argument("[Error]: Invalid local node index");
    }

    template<class ScalarType>
    typename Rectangle4P<ScalarType>::VectorType Rectangle4P<ScalarType>::toLocalPos(VectorType globalPos) const {
        return divide(ScalarType(2) * globalPos - bottomLeft - topRight, topRight - bottomLeft);
    }

    template<class ScalarType>
    typename Rectangle4P<ScalarType>::VectorType Rectangle4P<ScalarType>::toGlobalPos(VectorType localPos) const {
        return (bottomLeft + topRight + hadamard(bottomLeft - topRight, localPos)) * ScalarType(0.5);
    }

    template<class ScalarType>
    void Rectangle4P<ScalarType>::swap(Rectangle4P& elem) noexcept {
        Base::swap_base(elem);
        bottomLeft.swap(elem.bottomLeft);
        topRight.swap(elem.topRight);
    }

    template<class ScalarType>
    ScalarType Rectangle4P<ScalarType>::baseFunc(size_t localNode, VectorType p) {
        switch (localNode) {
            case BottomLeft:
                return (ScalarType(1) - p[0]) * (ScalarType(1) - p[1]) * ScalarType(0.25);
            case BottomRight:
                return (ScalarType(1) + p[0]) * (ScalarType(1) - p[1]) * ScalarType(0.25);
            case TopLeft:
                return (ScalarType(1) - p[0]) * (ScalarType(1) + p[1]) * ScalarType(0.25);
            case TopRight:
                return (ScalarType(1) + p[0]) * (ScalarType(1) + p[1]) * ScalarType(0.25);
        }
        throw std::invalid_argument("[Error]: Invalid local node index");
    }

    template<class ScalarType>
    ScalarType Rectangle4P<ScalarType>::dBase_dr(size_t localNode, [[maybe_unused]] VectorType p){
        switch (localNode) {
            case BottomLeft:
            case TopLeft:
                return -0.25;
            case BottomRight:
            case TopRight:
                return 0.25; 
        }
        throw std::invalid_argument("[Error]: Invalid local node index");
    }

    template<class ScalarType>
    ScalarType Rectangle4P<ScalarType>::dBase_ds(size_t localNode, [[maybe_unused]] VectorType p) {
        switch (localNode) {
            case BottomLeft:
            case BottomRight:
                return -0.25;
            case TopLeft:
            case TopRight:
                return 0.25; 
        }
        throw std::invalid_argument("[Error]: Invalid local node index");
    }

    template<class ScalarType>
    typename Rectangle4P<ScalarType>::VectorType Rectangle4P<ScalarType>::grad(size_t localNode, VectorType p) {
        return {dBase_dr(localNode, p), dBase_ds(localNode, p)};
    }
}