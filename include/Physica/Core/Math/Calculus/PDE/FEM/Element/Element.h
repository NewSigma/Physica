/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

namespace Physica {
    namespace Internal {
        template<class ElementType, class Functor, int DeltaOrder>
        struct GaussIntegral;
    }

    template<class Derived>
    class Element : public CRTP<Element<Derived>> {
        using This = Element<Derived>;
        using Base = CRTP<This>;
    public:
        using ScalarType = Traits<Derived>::ScalarType;
        constexpr static unsigned int Dim = Traits<Derived>::Dim;
        constexpr static unsigned int Order = Traits<Derived>::Order;
        constexpr static unsigned int DegreeOfFreedom = Traits<Derived>::DegreeOfFreedom;
        using IndexArray = Array<size_t, DegreeOfFreedom>;
        using VectorType = DenseVector<ScalarType, Dim>;
        using MatrixType = Traits<Derived>::MatrixType;
        static_assert(!ScalarType::isComplex(), "[Error]: Complex scalar is not allowed");
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
        template<int DeltaOrder = 0>
        [[nodiscard]] static ScalarType gauss_integral(std::invocable<VectorType> auto fn);
    protected:
        Element() = default;
        Element(IndexArray globalNodes_);
        void swap(Element& __restrict elem) noexcept;
    };

    template<class Derived>
    Element<Derived>::Element(IndexArray globalNodes_) : globalNodes(std::move(globalNodes_)) {}

    template<class Derived>
    template<int DeltaOrder>
    auto Element<Derived>::gauss_integral(std::invocable<VectorType> auto fn) -> ScalarType {
        return Internal::GaussIntegral<Derived, decltype(fn), DeltaOrder>::run(fn);
    }

    template<class Derived>
    void Element<Derived>::swap(Element& __restrict elem) noexcept {
        assert(this != &elem && "[Error]: Self swap is likely a bug");
        globalNodes.swap(elem.globalNodes);
    }
}

namespace Physica {
    template<class T>
    class Traits<Element<T>> {
    public:
        using Derived = T;
    };
}
