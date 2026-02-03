/*
 * Copyright 2025-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/SparseLU.cuh"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.cuh"
#include "Test.h"

using namespace Physica;
using T = float32;

namespace {
    void testSolver() {
        using VectorType = Vector4D<T>;
        using MatrixType = DenseMatrix<T, MatrixMajor::Row, 4, 4>;
        const MatrixType A{
                {-0.000696013585639699,     0.816492585748236, 0.0216969440126965, -0.0884307621566726},
                {    0.691809621910274, -0.000696013585639699,  0.131671000379563, -0.0701048797366553},
                {  -0.0701048797366553,   -0.0884307621566726, -0.131016640264434,   0.788769710999288},
                {    0.131671000379563,    0.0216969440126965,  0.643819646681362,  -0.131016640264434}
        };
        const VectorType b{4.316511702487202E-1, 1.548712563601895E-2, 9.840637243791538E-1, 1.671684099146560E-1};
        const VectorType answer{0.06910464034803039, 0.6682416388355244, 0.5106380624075890, 1.413471488683768};
        device_obj<SparseLU<T>> lu((SparseMatrix<T>(A)));
        auto d_b = device_obj<MatrixND<T>>(b.toDeviceAsync());
        lu.solve(d_b);
        auto result = d_b.toHost();
        expect(vectorNear(result.col(0), answer, 1E-6));
    }

    void testLnAbsDet() {
        using MatrixType = DenseMatrix<T, MatrixMajor::Row, 8, 8>;
        const MatrixType A = MatrixType::random_normal<Random<MCG, 1234>>(8, 8);
        device_obj<SparseLU<T>> lu((SparseMatrix<T>(A)));
        expect(scalarNear(lu.lnAbsDet(), A.lnAbsDet(), 1E-5));
        expect(lu.sgndet() == A.sgndet());
    }
}

int main() {
    testSolver();
    testLnAbsDet();
    return 0;
}
