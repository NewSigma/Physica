/*
 * Copyright 2022-2024 Weibo He.
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

#include "Physica/Core/Math/Calculus/PDE/FEM/Mesh.h"
#include "Element.h"

namespace Physica::Core {
    template<Scalar T>
    class Triangle1 : public Element<Triangle1<T>> {
    public:
        using Base = Element<Triangle1<T>>;
        using typename Base::VectorType;
        using typename Base::MatrixType;
        using typename Base::IndexArray;
        using PosArray = Array<VectorType, Traits<Triangle1<T>>::NumPoint>;
    private:
        PosArray pos;
    public:
        Triangle1() = default;
        Triangle1(PosArray pos_, IndexArray globalNodes);
        ~Triangle1() = default;
        /* Operators */
        Triangle1& operator=(Triangle1 elem) noexcept;
        /* Operations */
        void swap(Triangle1& __restrict elem) noexcept;
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
        [[nodiscard]] static Mesh<Triangle1> rectangle(VectorType bottomLeft,
                                                       VectorType topRight,
                                                       size_t numSeparateX,
                                                       size_t numSeparateY);
    };
}

namespace Physica {
    template<Scalar T>
    class Traits<Core::Triangle1<T>> {
    public:
        constexpr static unsigned int Dim = 2;
        constexpr static unsigned int Order = 1;
        constexpr static unsigned int NumPoint = 3;
        constexpr static unsigned int DegreeOfFreedom = NumPoint * Order;
        using ScalarType = T;
        using MatrixType = Core::DenseMatrix<T, MatrixOption::Col | MatrixOption::Element, Dim, Dim>;
    };
}

#include "Triangle1Impl.h"
