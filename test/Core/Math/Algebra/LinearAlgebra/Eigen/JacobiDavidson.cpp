/*
 * Copyright 2022-2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseHermiteMatrix.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/TransIsingMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<MT19937>;

namespace {
    void test1() {
        using MatrixType = DenseMatrix<cfloat64, MatrixMajor::Row>;
        const MatrixType data = MatrixType::random_uniform<RandomSource>(64);
        const DenseHermiteMatrix<cfloat64> hermite = data + data.hermite();
        RandomSource::getInstance().reseed(std::mt19937::default_seed);

        EigenSolver<cfloat64> eig(hermite, false);
        eig.sort();
        JacobiDavidson<cfloat64> jd(hermite.getRow(), 48);
        jd.compute(hermite, VectorND<cfloat64>::random_uniform<RandomSource>(data.getRow()));
        jd.sort();

        expect<RandomSource>(vectorNear(jd.getEigenvalues(), eig.getEigenvalues().head(jd.getNumRequired()), 1E-12));
    }

    void test2() {
        constexpr int Dim = 1;
        constexpr int NumSite = 8;
        constexpr int NumLevel = 6;
        constexpr BoundaryCond BC = BoundaryCond::PBC;
        using T = float64;
        using Hamiltonian = TransIsingMatrix<float64, Dim, NumSite, BC>;
        RandomSource::getInstance().reseed(10000);

        auto Js = VectorND<T>::linspace(0, 1, 200);
        SquareLattice<Dim, BC> lattice({NumSite}, 1);
        Hamiltonian hamilton(T(1) - Js[1], Js[1], lattice);
        JacobiDavidson<T> jd(hamilton.getRow(), NumLevel);
        jd.compute(hamilton, VectorND<T>::template random_uniform<RandomSource>(hamilton.getRow()));
        jd.sort();

        EigenSolver<T> eig(hamilton, false);
        eig.sort();
        expect<RandomSource>(vectorNear(eig.getEigenvalues().reals().head(4), jd.getEigenvalues().head(4), 1E-12));
    }
}

int main() {
    test1();
    test2();
    return 0;
}
