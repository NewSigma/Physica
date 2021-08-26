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
#include "TestHelper.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Hessenburg.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"

using namespace Physica::Core;

int main() {
    {
        using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 4, 4>;
        const MatrixType mat{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        const MatrixType H_answer{{1, -5.38516, 0, 0},
                                {-16.5269, 33.8276, -6.75721, 0},
                                {-1.36458, 0.591256, -0.827586, 0},
                                {0, 0, 0, 0}};
        Hessenburg hess(mat);
        MatrixType H = hess.getMatrixH();
        if (!matrixNear(H, H_answer, 1E-5))
            return 1;
        MatrixType Q = hess.getMatrixQ();
        MatrixType A = (Q * H).compute() * Q.transpose();
        if (!matrixNear(A, mat, 1E-15))
            return 1;
    }
    {
        using MatrixType = DenseMatrix<Scalar<Double, false>, DenseMatrixOption::Column | DenseMatrixOption::Vector, 8, 8>;
        const MatrixType mat{{2.5, 0, 0.866025, 0, 0.481812, 0, 0.318105, 0},
                            {0, 10.5, 0, 4.97494, 0, 3.06186, 0, 2.12934},
                            {0.866025, 0, 25.5, 0, 14.1869, 0, 9.36654, 0},
                            {0, 4.97494, 0, 49.5, 0, 30.4651, 0, 21.1867},
                            {0.481812, 0, 14.1869, 0, 84.5, 0, 55.789, 0},
                            {0, 3.06186, 0, 30.4651, 0, 132.5, 0, 92.1458},
                            {0.318105, 0, 9.36654, 0, 55.789, 0, 195.5, 0},
                            {0, 2.12934, 0, 21.1867, 0, 92.1458, 0, 275.5}};
        const MatrixType answer{{{2.5, -1.04083, 0, 0, 0, 0, 0, 0},
                                {-1.04083, 85.5001, -83.4666, 0, 0, 0, 0, 0},
                                {1.94772E-17, -83.4666, 161.216, 32.5875, 0, 0, 0, 0},
                                {2.51965E-17, -1.21897E-15, 32.5875, 58.7842, 4.32912E-14, 0, 0, 0},
                                {-2.45499E-16, 1.6109E-14, 5.8311E-15, 3.99302E-14, 33.304, 70.8353, 0, 0},
                                {7.73537E-17, -5.07573E-15, -1.41081E-15, 2.65133E-15, 70.8353, 297.412, 51.6424, 0},
                                {-7.67151E-17, 5.03383E-15, -3.80293E-15, 7.07156E-15, 1.05757E-14, 51.6424, 52.6065, 21.9094},
                                {2.83378E-17, -1.85945E-15, -1.80824E-16, 1.41259E-15, -7.22757E-15, 2.84217E-14, 21.9094, 84.6774}}};
        Hessenburg hess(mat);
        MatrixType result = hess.getMatrixH();
        if (!matrixNear(result, answer, 1E-5))
            return 1;

    }
    return 0;
}
