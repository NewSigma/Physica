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

namespace Physica::Core {
    template<class T>
    Mesh<T>::Mesh(size_t numElement, size_t numNode)
            : elements(numElement)
            , coeffs(numNode)
            , nodeTypes(numNode, NodeType::Free) {}

    template<class T>
    Mesh<T>& Mesh<T>::operator=(Mesh mesh) noexcept {
        this->swap(mesh);
        return *this;
    }

    template<class T>
    typename Mesh<T>::ScalarType Mesh<T>::operator()(VectorType p) const {
        for (const auto& elem : elements) {
            if (elem.contains(p)) {
                const auto& localNodes = elem.getNodes();
                const VectorType localPos = elem.toLocalPos(p);
                ScalarType result = 0;
                for (size_t i = 0; i < ElementType::getNumNodes(); ++i) {
                    const size_t globalNode = localNodes[i];
                    result += coeffs[globalNode] * ElementType::baseFunc(i, localPos);
                }
                return result;
            }
        }
        throw std::invalid_argument("[Error]: Accessing point outside domain of definition");
    }
    /**
     * \tparam Detector
     * A functor defined as
     * bool Detector(VectorType)
     * \returns true if the given node is on the boundary
     * 
     * \tparam Conditioner
     * A functor defined as
     * ScalarType Conditioner(VectorType)
     * \returns the dirichlet boundary value
     */
    template<class T>
    template<class Detector, class Conditioner>
    void Mesh<T>::addDirichletBoundary(Detector detector, Conditioner conditioner) {
        for (const auto& elem : elements) {
            for (size_t i = 0; i < ElementType::getNumNodes(); ++i) {
                const VectorType pos = elem.getNodePos(i);
                const bool isOnBound = detector(pos);
                if (isOnBound) {
                    const size_t node = elem.getNodes()[i];
                    coeffs[node] = conditioner(pos);
                    nodeTypes[node] = NodeType::Dirichlet;
                }
            }
        }
    }

    template<class T>
    size_t Mesh<T>::getNumFreeNodes() const {
        size_t num = 0;
        for (auto type : nodeTypes)
            num += type == NodeType::Free;
        return num;
    }

    template<class T>
    void Mesh<T>::setElem(ElementType elem, size_t index) {
        elements[index].swap(elem);
    }

    template<class T>
    Utils::Array<typename Mesh<T>::VectorType> Mesh<T>::getNodes() const {
        Utils::Array<typename Mesh<T>::VectorType> result(getNumNodes());
        for (const auto& element : elements) {
            size_t local_index = 0;
            for (size_t node : element.getNodes()) {
                result[node] = element.getNodePos(local_index);
                ++local_index;
            }
        }
        return result;
    }

    template<class T>
    void Mesh<T>::swap(Mesh& mesh) noexcept {
        elements.swap(mesh.elements);
        coeffs.swap(mesh.coeffs);
        nodeTypes.swap(nodeTypes);
    }

    template<class ScalarType>
    Mesh<Rectangle4P<ScalarType>> rectangle(Vector<ScalarType, 2> bottomLeft,
                                            Vector<ScalarType, 2> topRight,
                                            size_t numElementX,
                                            size_t numElementY) {
        using VectorType = Vector<ScalarType, 2>;
        const size_t numNodeX = numElementX + 1;
        const size_t numNodeY = numElementY + 1;
        Mesh<Rectangle4P<ScalarType>> mesh(numElementX * numElementY, numNodeX * numNodeY);
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
                mesh.setElem(Rectangle4P<ScalarType>(p, p + diagnal, nodeBottomLeft, nodeBottomRight, nodeTopLeft, nodeTopRight),
                             nextElem++);
                p[0] += xPerElem;
            }
            p[1] += yPerElem;
            p[0] = bottomLeft[0];
        }
        return mesh;
    }
}
