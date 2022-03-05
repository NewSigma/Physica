/*
 * Copyright 2021 WeiBo He.
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
#include "Physica/Utils/TestHelper.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearEquations.h"

using namespace Physica::Utils;
using namespace Physica::Core;

int main() {
    using ScalarType = Scalar<Double, false>;
    using MatrixType = DenseMatrix<ScalarType, DenseMatrixOption::Column | DenseMatrixOption::Row, 4, 5>;
    using VectorType = Vector<ScalarType, 4>;
    MatrixType A{{0.691809621910274, -0.000696013585639699, 0.131671000379563, -0.0701048797366553, 4.316511702487202E-1},
                {-0.000696013585639699, 0.816492585748236, 0.0216969440126965, -0.0884307621566726, 1.548712563601895E-2},
                {0.131671000379563, 0.0216969440126965, 0.643819646681362, -0.131016640264434, 9.840637243791538E-1},
                {-0.0701048797366553, -0.0884307621566726, -0.131016640264434, 0.788769710999288, 1.671684099146560E-1}};
    VectorType answer{0.379822910240522, 0.0329724647167976, 1.55292884169613, 0.507335846409689};
    {
        LinearEquations equ(A);
        equ.gaussJordanPartial();
        if (!vectorNear(equ.getSolution(), answer, 1E-14))
            return 1;
    }
    {
        LinearEquations equ(A);
        equ.gaussJordanComplete();
        if (!vectorNear(equ.getSolution(), answer, 1E-13))
            return 1;
    }
    {
        LinearEquations equ(A);
        equ.gaussEliminationPartial();
        if (!vectorNear(equ.getSolution(), answer, 1E-13))
            return 1;
    }
    {
        LinearEquations equ(A);
        equ.gaussEliminationComplete();
        if (!vectorNear(equ.getSolution(), answer, 1E-14))
            return 1;
    }
    {
        LinearEquations equ(A);
        equ.lu();
        if (!vectorNear(equ.getSolution(), answer, 1E-14))
            return 1;
    }
    {
        LinearEquations equ(A);
        equ.conjugateGradient();
        if (!vectorNear(equ.getSolution(), answer, 1E-14))
            return 1;
    }
    return 0;
}
