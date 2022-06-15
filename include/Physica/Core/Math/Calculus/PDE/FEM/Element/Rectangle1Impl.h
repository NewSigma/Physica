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
    Rectangle1<ScalarType>::Rectangle1(VectorType bottomLeft_,
                                       VectorType topRight_,
                                       IndexArray globalNodes)
            : Base(std::move(globalNodes)), bottomLeft(bottomLeft_), topRight(topRight_) {}

    template<class ScalarType>
    Rectangle1<ScalarType>& Rectangle1<ScalarType>::operator=(Rectangle1<ScalarType> elem) noexcept {
        swap(elem);
        return *this;
    }

    template<class ScalarType>
    typename Rectangle1<ScalarType>::MatrixType Rectangle1<ScalarType>::jacobi([[maybe_unused]] VectorType localPos) const {
        return MatrixType{(topRight[0] - bottomLeft[0]) * 0.5, 0, 0, (topRight[1] - bottomLeft[1]) * 0.5};
    }

    template<class ScalarType>
    typename Rectangle1<ScalarType>::MatrixType Rectangle1<ScalarType>::inv_jacobi([[maybe_unused]] VectorType globalPos) const {
        return MatrixType{ScalarType::Two() / (topRight[0] - bottomLeft[0]), 0, 0, ScalarType::Two() / (topRight[1] - bottomLeft[1])};
    }

    template<class ScalarType>
    bool Rectangle1<ScalarType>::contains(const VectorType& point) const {
        return bottomLeft[0] <= point[0]
            && point[0] <= topRight[0]
            && bottomLeft[1] <= point[1]
            && point[1] <= topRight[1];
    }

    template<class ScalarType>
    typename Rectangle1<ScalarType>::VectorType Rectangle1<ScalarType>::getNodePos(size_t localNode) const {
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
    typename Rectangle1<ScalarType>::VectorType Rectangle1<ScalarType>::toLocalPos(VectorType globalPos) const {
        return divide(ScalarType(2) * globalPos - bottomLeft - topRight, topRight - bottomLeft);
    }

    template<class ScalarType>
    typename Rectangle1<ScalarType>::VectorType Rectangle1<ScalarType>::toGlobalPos(VectorType localPos) const {
        return (bottomLeft + topRight + hadamard(bottomLeft - topRight, localPos)) * ScalarType(0.5);
    }

    template<class ScalarType>
    void Rectangle1<ScalarType>::swap(Rectangle1& elem) noexcept {
        Base::swap(elem);
        bottomLeft.swap(elem.bottomLeft);
        topRight.swap(elem.topRight);
    }

    template<class ScalarType>
    ScalarType Rectangle1<ScalarType>::baseFunc(size_t localNode, VectorType p) {
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
    ScalarType Rectangle1<ScalarType>::dBase_dr(size_t localNode, [[maybe_unused]] VectorType p){
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
    ScalarType Rectangle1<ScalarType>::dBase_ds(size_t localNode, [[maybe_unused]] VectorType p) {
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
    typename Rectangle1<ScalarType>::VectorType Rectangle1<ScalarType>::grad(size_t localNode, VectorType p) {
        return {dBase_dr(localNode, p), dBase_ds(localNode, p)};
    }

    template<class ScalarType>
    template<class Functor>
    ScalarType Rectangle1<ScalarType>::integral(Functor func) {
        constexpr static unsigned int Factor = 1U << Base::Dim;
        return func(VectorType{0, 0}) * Factor;
    }

    template<class ScalarType>
    Mesh<Rectangle1<ScalarType>> Rectangle1<ScalarType>::rectangle(VectorType bottomLeft,
                                                                   VectorType topRight,
                                                                   size_t numElementX,
                                                                   size_t numElementY) {
        using VectorType = Vector<ScalarType, 2>;
        const size_t numNodeX = numElementX + 1;
        const size_t numNodeY = numElementY + 1;
        Mesh<Rectangle1<ScalarType>> mesh(numElementX * numElementY, numNodeX * numNodeY);
        const ScalarType xPerElem = (topRight[0] - bottomLeft[0]) / ScalarType(numElementX);
        const ScalarType yPerElem = (topRight[1] - bottomLeft[1]) / ScalarType(numElementY);
        const VectorType diagnal{xPerElem, yPerElem};

        VectorType p = bottomLeft;
        size_t nextElem = 0;
        for (size_t y = 0; y < numElementY; ++y) {
            for (size_t x = 0; x < numElementX; ++x) {
                size_t nodeBottomLeft = numNodeX * y + x;
                size_t nodeBottomRight = nodeBottomLeft + 1;
                size_t nodeTopLeft = nodeBottomLeft + numNodeX;
                size_t nodeTopRight = nodeBottomRight + numNodeX;
                mesh.setElem(Rectangle1<ScalarType>(p, p + diagnal, {nodeBottomLeft, nodeBottomRight, nodeTopRight, nodeTopLeft}),
                             nextElem++);
                p[0] += xPerElem;
            }
            p[1] += yPerElem;
            p[0] = bottomLeft[0];
        }
        return mesh;
    }
}
