/*
 * Copyright 2022-2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/Diff.h"

using namespace Physica;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

namespace {
    template<Matrix M, Backend B>
    void eigenTestImpl(const M& mat, double precision) {
        using ScalarType = M::ScalarType;
        using VectorType = DenseVector<ScalarType, M::RowAtCompile>;
        using EigenvectorMatrix = SymmEigenSolver<ScalarType>::EigenvectorMatrix;

        auto solver = SymmEigenSolver<ScalarType>(mat.getRow(), true);
        if constexpr (B == Backend::MKL)
            solver.compute_mkl(mat);
        else
            solver.compute_base(mat);
        solver.sort();

        const size_t order = mat.getRow();
        EigenvectorMatrix eigenvectors = solver.getEigenvectors();
        for (size_t i = 0; i < order; ++i) {
            if (i > 1 && solver.getEigenvalues()[i - 1] > solver.getEigenvalues()[i])
                exit(EXIT_FAILURE);

            VectorType v1 = mat * eigenvectors.col(i);
            VectorType v2 = eigenvectors.col(i) * ScalarType(solver.getEigenvalues()[i]);
            if (!vectorNear(v1, v2, precision))
                exit(EXIT_FAILURE);
        }
    }

    void eigenTest(const Matrix auto& m, double prec) {
        using M = std::remove_cvref_t<decltype(m)>;
        eigenTestImpl<M, Backend::Base>(m, prec);

        if constexpr (HasMKL() && !Diffable<M>)
            eigenTestImpl<M, Backend::MKL>(m, prec);
    }
}

int main() {
    {
        using MatrixType = DenseMatrix<float64, MatrixOption::Col, 3, 3>;
        const MatrixType data{
                {-0.590316, -2.195140, -2.374630},
                {-1.250060, -0.297493,  1.403490},
                { 0.517063, -0.956614, -0.920775}
        };
        const MatrixType mat = data + data.transpose();
        eigenTest(mat, 1E-14);

        using T = Diff<float64, DiffMode::Forward, 1>;
        using DiffMatrix = DenseMatrix<T, MatrixOption::Col, 3, 3>;
        const DiffMatrix mat1 = mat;
        eigenTest(mat1, 1E-14);
    }
    {
        using MatrixType = DenseSymmMatrix<float64>;
        const auto mat = MatrixType::random_uniform<RandomSource>(8);
        eigenTest(mat, 1E-12);
    }
}
