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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica::Core {
    template<class Derived>
    class Element : public Utils::CRTPBase<Derived> {
        using Base = Utils::CRTPBase<Derived>;
    public:
        using ScalarType = typename Internal::Traits<Derived>::ScalarType;
        constexpr static unsigned int Dim = Internal::Traits<Derived>::Dim;
        constexpr static unsigned int Order = Internal::Traits<Derived>::Order;
        constexpr static unsigned int DegreeOfFreedom = Internal::Traits<Derived>::DegreeOfFreedom;
        using IndexArray = Utils::Array<size_t, DegreeOfFreedom>;
        using VectorType = Vector<ScalarType, Dim>;
        using MatrixType = typename Internal::Traits<Derived>::MatrixType;
    protected:
        IndexArray globalNodes;
    public:
        /* Getters */
        [[nodiscard]] MatrixType jacobi(VectorType localPos) const { return Base::getDerived().jacobi(localPos); }
        [[nodiscard]] MatrixType inv_jacobi(VectorType globalPos) const { return Base::getDerived().inv_jacobi(globalPos); }
        [[nodiscard]] bool contains(const VectorType& point) const { return Base::getDerived().contains(point); }
        [[nodiscard]] const IndexArray& getGlobalNodes() const { return globalNodes; }
        [[nodiscard]] VectorType getNodePos(size_t localNode) const { return Base::getDerived().getNodePos(localNode); }
        [[nodiscard]] VectorType toLocalPos(VectorType globalPos) const { return Base::getDerived().toLocalPos(globalPos); }
        [[nodiscard]] VectorType toGlobalPos(VectorType localPos) const { return Base::getDerived().toLocalPos(localPos); }
        [[nodiscard]] constexpr static size_t getNumNodes() { return DegreeOfFreedom; }
    protected:
        Element() = default;
        Element(IndexArray globalNodes_);
        void swap(Element& elem) noexcept;
    };

    template<class Derived>
    Element<Derived>::Element(IndexArray globalNodes_) : globalNodes(std::move(globalNodes_)) {}

    template<class Derived>
    void Element<Derived>::swap(Element& elem) noexcept {
        globalNodes.swap(elem.globalNodes);
    }
}
