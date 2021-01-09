/*
 * Copyright 2020 WeiBo He.
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
        //typedef Matrix<Scalar, MatrixType::ElementColumn> Matrix;
    public:
        ~Element();
        /* Operations */
        virtual Scalar determinate(const Scalar& s1, const Scalar& s2) = 0;
        virtual Scalar shapePartialS1(size_t shapeIndex, const Scalar& s1, const Scalar& s2) = 0;
        virtual Scalar shapePartialS2(size_t shapeIndex, const Scalar& s1, const Scalar& s2) = 0;
        //virtual Matrix jacobi(const Scalar& s1, const Scalar& s2) = 0;
        //virtual Matrix inverseJacobi(const Scalar& s1, const Scalar& s2) = 0;
    protected:
        explicit Element(size_t nodesCount);
    };
}

#endif
