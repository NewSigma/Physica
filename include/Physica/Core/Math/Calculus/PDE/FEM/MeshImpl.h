/*
 * Copyright 2022-2025 Weibo He.
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

#include "Mesh.h"

namespace Physica {
    template<class T>
    Mesh<T>::Mesh(size_t numElement, size_t numNode)
            : elements(numElement)
            , coeffs(numNode)
            , nodeTypes(numNode, NodeType::Free) {}

    template<class T>
    auto Mesh<T>::operator()(VectorType p) const -> ScalarType {
        for (const auto& elem : elements) {
            if (elem.contains(p)) {
                const auto& globalNodes = elem.getGlobalNodes();
                const VectorType localPos = elem.toLocalPos(p);
                ScalarType result = 0;
                for (size_t i = 0; i < ElementType::getNumNodes(); ++i) {
                    const size_t globalNode = globalNodes[i];
                    result += coeffs[globalNode] * ElementType::baseFunc(i, localPos);
                }
                return result;
            }
        }
        throw std::invalid_argument("[Error]: Accessing point outside domain of definition");
    }

    template<class T>
    void Mesh<T>::addDirichletBoundary(std::invocable<VectorType> auto detector, std::invocable<VectorType> auto conditioner) {
        for (const auto& elem : elements) {
            for (size_t i = 0; i < ElementType::getNumNodes(); ++i) {
                const VectorType pos = elem.getNodePos(i);
                const bool isOnBound = detector(pos);
                if (isOnBound) {
                    const size_t node = elem.getGlobalNodes()[i];
                    coeffs[node] = conditioner(pos);
                    nodeTypes[node] = NodeType::Dirichlet;
                }
            }
        }
    }

    template<class T>
    void Mesh<T>::swap(Mesh& __restrict mesh) noexcept {
        assert(this != &mesh && "[Error]: Self swap is likely a bug");
        elements.swap(mesh.elements);
        coeffs.swap(mesh.coeffs);
        nodeTypes.swap(nodeTypes);
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
    auto Mesh<T>::getNodes() const -> Array<VectorType> {
        Array<VectorType> result(getNumNodes());
        for (const auto& element : elements) {
            size_t local_index = 0;
            for (size_t node : element.getGlobalNodes()) {
                result[node] = element.getNodePos(local_index);
                ++local_index;
            }
        }
        return result;
    }
}
