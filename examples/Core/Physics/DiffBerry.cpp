/*
 * Copyright 2025 Weibo He.
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
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/ForwardEigenSolver.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/PauliMatrix.h"
#include "Physica/Core/Math/Calculus/Integrate/Vegas.h"

using namespace Physica;
using T = float64;
using dfloat = Diff<T, DiffMode::Forward>;
using dcfloat = Diff<cfloat64, DiffMode::Forward>;

namespace {
    T curv(Vector2D<T> k0, int band) {
        auto diag = [k0, band](int var) -> Vector2D<cfloat64> {
            Vector2D<dfloat> k = k0;
            k[var].grad() = 1;
            DenseMatrix<dcfloat> hamilton = sin(k[0]) * PauliMatrix<T, PauliIndex::X>{}
                                          + T(3) * sin(k[1]) * PauliMatrix<T, PauliIndex::Y>{}
                                          + (T(1) - cos(k[0]) - cos(k[1])) * PauliMatrix<T, PauliIndex::Z>{};
            EigenSolver<dcfloat, 2> eig(hamilton, true);
            eig.sort();
            return eig.getEigenvectors().col(band).grads();
        };
        return T(-2) *(diag(0).conjugate() * diag(1)).imag();
    }

    void chern(int band) {
        auto vegas = Vegas<T>({-std::numbers::pi, -std::numbers::pi}, {std::numbers::pi, std::numbers::pi}, 100, 10000);
        vegas.integral<Random<>>([band](Vector2D<T> k0) {
            return curv(k0, band);
        });
        T factor = reciprocal(MathConst<T>::pi * 2);
        std::cout << std::format("{}({})\n", vegas.calcMean() * factor, vegas.calcDevia() * factor);
    }
}

int main() {
    chern(0);
    chern(1);
    return 0;
}
