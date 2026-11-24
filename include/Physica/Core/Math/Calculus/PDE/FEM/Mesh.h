/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Geometry/Point.h"

namespace Physica {
    enum class NodeType : char {
        Free,
        Dirichlet
    };

    template<class T>
    class Mesh {
    public:
        using ElementType = T;
        using ScalarType = ElementType::ScalarType;
        using VectorType = ElementType::VectorType;
    private:
        Array<ElementType> elements;
        VectorND<ScalarType> coeffs;
        Array<NodeType> nodeTypes;
    public:
        Mesh(size_t numElement, size_t numNode);
        Mesh(const Mesh&) = default;
        Mesh(Mesh&&) noexcept = default;
        ~Mesh() = default;
        /* Operators */
        Mesh& operator=(Mesh obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] ScalarType operator()(VectorType p) const;
        /* Operations */
        void addDirichletBoundary(std::invocable<VectorType> auto detector, std::invocable<VectorType> auto conditioner);
        void swap(Mesh& __restrict mesh) noexcept;
        /* Getters */
        [[nodiscard]] size_t getNumElems() const { return elements.getLength(); }
        [[nodiscard]] size_t getNumNodes() const { return coeffs.getLength(); }
        [[nodiscard]] const Array<ElementType>& getElements() const { return elements; }
        [[nodiscard]] VectorND<ScalarType>& getCoeffs() { return coeffs; }
        [[nodiscard]] const VectorND<ScalarType>& getCoeffs() const { return coeffs; }
        [[nodiscard]] const Array<NodeType>& getNodeTypes() const { return nodeTypes; }
        [[nodiscard]] size_t getNumFreeNodes() const;
        [[nodiscard]] Array<VectorType> getNodes() const;
        /* Setters */
        void setElem(ElementType elem, size_t index);
    };

    template<class T>
    void swap(Mesh<T>& mesh1, Mesh<T>& mesh2) noexcept {
        mesh1.swap(mesh2);
    }
}

#include "MeshImpl.h"
