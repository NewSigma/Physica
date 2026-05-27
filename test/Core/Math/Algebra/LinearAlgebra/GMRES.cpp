/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearSystem/GMRES.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using Tc = cfloat64;

namespace {
    void simple() noexcept {
        /* Real */ {
            DenseMatrix<T, MatrixMajor::Row, 3, 3> A{
                {3,  2, 0},
                {1, -1, 0},
                {0,  5, 1}
            };
            const VectorND<T> b{2, 4, -1};
            auto x = b;
            auto gmres = GMRES<T>(3, 3);
            gmres.solve(A, x);
            expect(vectorNear(A * x, b, 1E-15));
        }
        /* Complex */ {
            DenseMatrix<Tc, MatrixMajor::Row, 3, 3> A{
                {{2, 1},      -1,  0},
                {    -1, {2, -1}, -1},
                {     0,  {1, 2},  3}
            };
            const VectorND<Tc> b{{-1, 4}, {-1, -5}, {13, 3}};
            auto x = b;
            auto gmres = GMRES<Tc>(3, 3);
            gmres.solve(A, x);
            expect(vectorNear(A * x, b, 1E-14));
        }
        /* FOM fails */ {
            DenseMatrix<T, MatrixMajor::Row, 3, 3> A{
                {2.5789473684210527,  2.1951219512195124, 0},
                {0.7894736842105263, -1.1951219512195121, 0},
                {                 0,   5.487804878048781, 1}
            };
            const VectorND<T> b{1, 1, 1};
            auto x = b;
            auto gmres = GMRES<T>(3, 3);
            gmres.solve(A, x);
            expect(vectorNear(A * x, b, 1E-15));
        }
    }

    void eigen() noexcept {
        /* Test that the RHS is an eigenvector */ {
            DenseMatrix<T, MatrixMajor::Row, 2, 2> A{
                {1, 2},
                {3, 0}
            };
            VectorND<T> x{1, 1};
            auto gmres = GMRES<T>(2, 2);
            gmres.setIterationLimit(1);
            gmres.solve(A, x);
            expect(vectorNear(x, VectorND<T>{1.0 / 3, 1.0 / 3}, 1E-15));
        }
        /* Test that the RHS is an eigenvector with eigenvalue 1 */ {
            DenseMatrix<T, MatrixMajor::Row, 3, 3> A{
                {3,  2, 0},
                {1, -1, 0},
                {0,  5, 1}
            };
            VectorND<T> x{0, 0, 1};
            auto gmres = GMRES<T>(3, 3);
            gmres.solve(A, x);
            expect(vectorNear(A * x, x, 1E-15));
        }
    }

    void exact() {
        // Test GMRES do not iterate if we provide exact solution
        DenseMatrix<T, MatrixMajor::Row, 3, 3> A{
            {1, 2, 0},
            {2, 1, 2},
            {0, 2, 1}
        };
        VectorND<T> b{3, 5, 3};
        VectorND<T> x{1, 1, 1};
        auto gmres = GMRES<T>(3, 3);
        gmres.solve(A, b, x);
        expect(gmres.getIteration() == 0);
    }
}

int main() {
    simple();
    eigen();
    exact();
    return 0;
}
