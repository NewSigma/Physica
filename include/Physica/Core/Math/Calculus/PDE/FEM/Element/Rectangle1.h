/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Calculus/PDE/FEM/Mesh.h"
#include "Element.h"

namespace Physica {
    template<Scalar T>
    class Rectangle1 : public Element<Rectangle1<T>> {
    public:
        using Base = Element<Rectangle1<T>>;
        using typename Base::VectorType;
        using typename Base::MatrixType;
        using typename Base::IndexArray;
    private:
        enum {
            BottomLeft = 0,
            BottomRight = 1,
            TopRight = 2,
            TopLeft = 3
        };

        VectorType bottomLeft;
        VectorType topRight;
    public:
        Rectangle1() = default;
        Rectangle1(VectorType bottomLeft_, VectorType topRight_, IndexArray globalNodes);
        Rectangle1(const Rectangle1&) = default;
        Rectangle1(Rectangle1&&) noexcept = default;
        ~Rectangle1() = default;
        /* Operators */
        Rectangle1& operator=(Rectangle1 elem) noexcept;
        /* Operations */
        void swap(Rectangle1& __restrict elem) noexcept;
        /* Getters */
        [[nodiscard]] MatrixType jacobi([[maybe_unused]] VectorType localPos) const;
        [[nodiscard]] MatrixType inv_jacobi([[maybe_unused]] VectorType globalPos) const;
        [[nodiscard]] bool contains(const VectorType& point) const;
        [[nodiscard]] VectorType getNodePos(size_t localNode) const;
        [[nodiscard]] VectorType toLocalPos(VectorType globalPos) const;
        [[nodiscard]] VectorType toGlobalPos(VectorType localPos) const;
        /* Static members */
        [[nodiscard]] static T baseFunc(size_t localNode, VectorType p);
        [[nodiscard]] static T dBase_dr(size_t localNode, [[maybe_unused]] VectorType p);
        [[nodiscard]] static T dBase_ds(size_t localNode, [[maybe_unused]] VectorType p);
        [[nodiscard]] static VectorType grad(size_t localNode, [[maybe_unused]] VectorType p);
        [[nodiscard]] static Mesh<Rectangle1<T>> rectangle(VectorType bottomLeft,
                                                                    VectorType topRight,
                                                                    size_t numElementX,
                                                                    size_t numElementY);
    };
}

namespace Physica {
    template<Scalar T>
    class Traits<Rectangle1<T>> {
    public:
        constexpr static unsigned int Dim = 2;
        constexpr static unsigned int Order = 1;
        constexpr static unsigned int NumPoint = 4;
        constexpr static unsigned int DegreeOfFreedom = NumPoint * Order;
        using ScalarType = T;
        using MatrixType = DenseMatrix<T, MatrixMajor::Col, Dim, Dim>;
    };
}

#include "Rectangle1Impl.h"
