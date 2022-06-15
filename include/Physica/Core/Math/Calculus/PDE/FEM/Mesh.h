/*
 * Copyright 2020-2022 WeiBo He.
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

#include "Physica/Core/Math/Geometry/Point.h"

namespace Physica::Core {
    enum class NodeType : char {
        Free,
        Dirichlet
    };

    template<class T>
    class Mesh {
    public:
        using ElementType = T;
        using ScalarType = typename ElementType::ScalarType;
        using VectorType = typename ElementType::VectorType;
    private:
        Utils::Array<ElementType> elements;
        Vector<ScalarType> coeffs;
        Utils::Array<NodeType> nodeTypes;
    public:
        Mesh(size_t numElement, size_t numNode);
        Mesh(const Mesh&) = default;
        Mesh(Mesh&&) noexcept = default;
        ~Mesh() = default;
        /* Operators */
        Mesh& operator=(Mesh mesh) noexcept;
        [[nodiscard]] ScalarType operator()(VectorType p) const;
        template<class Detector, class Conditioner>
        void addDirichletBoundary(Detector detector, Conditioner conditioner);
        /* Getters */
        [[nodiscard]] size_t getNumElems() const { return elements.getLength(); }
        [[nodiscard]] size_t getNumNodes() const { return coeffs.getLength(); }
        [[nodiscard]] const Utils::Array<ElementType>& getElements() const { return elements; }
        [[nodiscard]] Vector<ScalarType>& getCoeffs() { return coeffs; }
        [[nodiscard]] const Vector<ScalarType>& getCoeffs() const { return coeffs; }
        [[nodiscard]] const Utils::Array<NodeType>& getNodeTypes() const { return nodeTypes; }
        [[nodiscard]] size_t getNumFreeNodes() const;
        [[nodiscard]] Utils::Array<VectorType> getNodes() const;
        /* Setters */
        void setElem(ElementType elem, size_t index);
        /* Helpers */
        void swap(Mesh& mesh) noexcept;
    };

    template<class T>
    inline void swap(Mesh<T>& mesh1, Mesh<T>& mesh2) noexcept {
        mesh1.swap(mesh2);
    }
}

#include "MeshImpl.h"
