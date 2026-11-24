/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
