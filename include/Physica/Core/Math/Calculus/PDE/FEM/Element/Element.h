/*
 * Copyright 2020-2021 WeiBo He.
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
#ifndef PHYSICA_ELEMENT_H
#define PHYSICA_ELEMENT_H

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "AbstractElement.h"

namespace Physica::Core {
    template<size_t dim, class Scalar>
    class Element;

    template<class Scalar>
    class Element<2, Scalar> : public AbstractElement<2> {
    public:
        typedef DenseMatrix<Scalar, DenseMatrixType::ElementColumn> Matrix;
        /**
         * There are only 4 shape functions. This enum is designed for performance when using switch.
         */
        enum ShapeIndex {
            Shape1, Shape2, Shape3, Shape4
        }
    public:
        explicit Element(size_t nodesCount);
        ~Element();
        /* Operations */
        //Optimize: The following functions are frequantly used. Remove the virtual if possible.
        /**
         * \return
         * The determinate of jacobi matrix.
         * \param s1
         * Coordinate in bi-unit element.
         * \param s2
         * Coordinate in bi-unit element.
         */
        virtual Scalar determinate(const Poing<2>& p) final { return jacobi(p).inverse(); };
        /**
         * \return
         * The derivative of the \param shapeIndex.th shape function with respect to \param s1 at (s1, s2).
         */
        virtual Scalar shapePartialS1(size_t shapeIndex, const Poing<2>& p) = 0;
        /**
         * \return
         * The derivative of the \param shapeIndex.th shape function with respect to \param s2 at (s1, s2).
         */
        virtual Scalar shapePartialS2(size_t shapeIndex, const Poing<2>& p) = 0;
        virtual Matrix jacobi(const Poing<2>& p) = 0;
        virtual Matrix inverseJacobi(const Poing<2>& p) final { return jacobi(p).inverse(); }
    };
}

#endif
