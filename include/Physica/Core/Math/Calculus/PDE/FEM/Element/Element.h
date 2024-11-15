/*
 * Copyright 2020-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h>

namespace Physica::Core {
    namespace Internal {
        template<class ElementType, class Functor, int DeltaOrder>
        struct GaussIntegral;
    }

    template<class Derived>
    class Element : public CRTPBase<Element<Derived>> {
        using This = Element<Derived>;
        using Base = CRTPBase<This>;
    public:
        using ScalarType = typename Traits<Derived>::ScalarType;
        constexpr static unsigned int Dim = Traits<Derived>::Dim;
        constexpr static unsigned int Order = Traits<Derived>::Order;
        constexpr static unsigned int DegreeOfFreedom = Traits<Derived>::DegreeOfFreedom;
        using IndexArray = Array<size_t, DegreeOfFreedom>;
        using VectorType = DenseVector<ScalarType, Dim>;
        using MatrixType = typename Traits<Derived>::MatrixType;
        static_assert(!ScalarType::isComplex, "[Error]: Complex scalar is not allowed");
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
        /* Static members */
        template<class Functor, int DeltaOrder = 0>
        [[nodiscard]] static ScalarType gauss_integral(Functor func);
    protected:
        Element() = default;
        Element(IndexArray globalNodes_);
        void swap(Element& __restrict elem) noexcept;
    };

    template<class Derived>
    Element<Derived>::Element(IndexArray globalNodes_) : globalNodes(std::move(globalNodes_)) {}

    template<class Derived>
    template<class Functor, int DeltaOrder>
    typename Element<Derived>::ScalarType Element<Derived>::gauss_integral(Functor func) {
        return Internal::GaussIntegral<Derived, Functor, DeltaOrder>::run(func);
    }

    template<class Derived>
    void Element<Derived>::swap(Element& __restrict elem) noexcept {
        assert(this != &elem && "[Error]: Self swap is likely a bug");
        globalNodes.swap(elem.globalNodes);
    }
}

namespace Physica {
    template<class T>
    class Traits<Core::Element<T>> {
    public:
        using Derived = T;
    };
}
