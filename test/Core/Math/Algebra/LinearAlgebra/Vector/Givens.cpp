/*
 * Copyright 2021-2026 Weibo He.
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
#include "Physica/Core/Scalar/Complex.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;
constexpr float64 Prec = 1E-15;

namespace {
    template<Scalar T>
    void eliminateReal() {
        /* Single */ {
            const auto v = Vector2D<T>::template random_normal<RandomSource>(2);
            auto x = v;
            auto rotator = givens(x, 0, 1);
            applyGivensCol(rotator, x, 0, 1);
            expect<RandomSource>(x[1].value() < Prec);

            applyRowGivens(x, rotator, 0, 1);
            expect<RandomSource>(vectorNear(x, v, 5UL));
        }
        /* Row major */ {
            const auto m = DenseMatrix<T, MatrixMajor::Row>::template random_normal<RandomSource>(2, 16);
            auto x = m;
            auto rotator = givens(x.col(0), 0, 1);
            applyGivens(rotator, x, 0, 1);
            expect<RandomSource>(abs(x[1, 0].value()) < Prec);
        }
    }

    void eliminateComplex() {
        using Tc = cfloat64;
        /* Single */ {
            const auto v = Vector2D<Tc>::random_normal<RandomSource>(2);
            auto x = v;
            auto rotator = givens(x, 0, 1);
            applyGivensCol(rotator, x, 0, 1);
            expect<RandomSource>(x[1].norm() < Prec);

            applyRowGivens(x, rotator.conjugate(), 0, 1);
            expect<RandomSource>(vectorNear(x, v, 1E-14));
        }
        /* Row major */ {
            auto x = DenseMatrix<Tc, MatrixMajor::Row>::template random_normal<RandomSource>(2, 16);
            auto rotator = givens(x.col(0), 0, 1);
            applyGivens(rotator, x, 0, 1);
            expect<RandomSource>(x[1, 0].norm() < Prec);
        }
    }

    void emptyTest() {
        Vector2D<float64> v{1, 0};
        auto g = givens(v, 0, 1);
        expect(g[0] == float64(1));
        expect(g[1] == float64(0));
    }
}

int main() {
    eliminateReal<float64>();
    eliminateReal<Diff<float64, DiffMode::Forward, 1>>();
    eliminateComplex();
    emptyTest();
    return 0;
}
