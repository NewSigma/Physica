/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <iostream>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/ForwardEigenSolver.h"
#include "Physica/Core/Math/Calculus/Integrate/Vegas.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/PauliMatrix.h"

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
        return T(-2) * (diag(0).conjugate() * diag(1)).imag();
    }

    void chern(int band) {
        using Cube = Vegas<T>::Cube;
        auto vegas = Vegas<T>(Cube{{-std::numbers::pi, -std::numbers::pi}, {std::numbers::pi, std::numbers::pi}}, 10000);
        const auto process = vegas.integral<Random<>>([band](Vector2D<T> k0) {
            return curv(k0, band);
        }, 100);
        T factor = reciprocal(MathConst<T>::pi * 2);
        std::cout << std::format("{}({})\n", process.calcMean() * factor, process.calcDevia() * factor);
    }
}

int main() {
    chern(0);
    chern(1);
    return 0;
}
