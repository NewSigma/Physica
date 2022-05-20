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

#include "Element.h"

namespace Physica::Core {
    template<class ScalarType> class Rectangle4P;

    namespace Internal {
        template<class T>
        class Traits<Rectangle4P<T>> {
        public:
            constexpr static unsigned int Dim = 2;
            constexpr static unsigned int NumPoint = 4;
            constexpr static unsigned int Order = 1;
            constexpr static unsigned int DegreeOfFreedom = NumPoint * Order;
            using ScalarType = T;
            using MatrixType = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Element, Dim, Dim>;
        };
    }

    template<class ScalarType>
    class Rectangle4P : public Element<Rectangle4P<ScalarType>> {
    public:
        using Base = Element<Rectangle4P<ScalarType>>;
        using typename Base::VectorType;
        using MatrixType = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Element, 2, 2>;
    private:
        enum {
            BottomLeft = 0,
            BottomRight = 1,
            TopLeft = 2,
            TopRight = 3
        };

        VectorType bottomLeft;
        VectorType topRight;
    public:
        Rectangle4P() = default;
        Rectangle4P(VectorType bottomLeft_,
                    VectorType topRight_,
                    size_t nodeBottomLeft,
                    size_t nodeBottomRight,
                    size_t nodeTopLeft,
                    size_t nodeTopRight);
        ~Rectangle4P() = default;
        /* Operators */
        Rectangle4P& operator=(Rectangle4P elem) noexcept;
        /* Getters */
        [[nodiscard]] MatrixType jacobi([[maybe_unused]] VectorType localPos) const;
        [[nodiscard]] MatrixType inv_jacobi([[maybe_unused]] VectorType globalPos) const;
        [[nodiscard]] bool contains(const VectorType& point) const;
        [[nodiscard]] VectorType getNodePos(size_t localNode) const;
        [[nodiscard]] VectorType toLocalPos(VectorType globalPos) const;
        [[nodiscard]] VectorType toGlobalPos(VectorType localPos) const;
        /* Helpers */
        void swap(Rectangle4P& elem) noexcept;
        /* Static members */
        [[nodiscard]] static ScalarType baseFunc(size_t localNode, VectorType p);
        [[nodiscard]] static ScalarType dBase_dr(size_t localNode, [[maybe_unused]] VectorType p);
        [[nodiscard]] static ScalarType dBase_ds(size_t localNode, [[maybe_unused]] VectorType p);
        [[nodiscard]] static VectorType grad(size_t localNode, [[maybe_unused]] VectorType p);
    };
}

#include "Rectangle4PImpl.h"
