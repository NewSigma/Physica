/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/Tridiagonalization.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Scalar/Complex.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<MT19937, std::mt19937::default_seed>;

namespace {
    template<Matrix M>
    void decomp(const M& source, double tolerance) {
        using ScalarType = M::ScalarType;
        using GeneralMatrix = DenseMatrix<ScalarType>;
        Tridiagonalization<ScalarType> tri(source);
        GeneralMatrix T = tri.getMatrixT();
        GeneralMatrix Q = tri.getMatrixQ();
        GeneralMatrix A = (Q * T).compute() * Q.hermite();
        expect(matrixNear(A, source, tolerance));
    }
}

int main() {
    using RealType = float64;
    {
        using MatrixType = DenseMatrix<RealType, MatrixMajor::Col, 4, 4>;
        const MatrixType temp{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};
        const MatrixType mat = temp + temp.transpose();
        decomp(mat, 1E-14);
    }
    {
        using MatrixType = DenseSymmMatrix<RealType>;
        const auto mat = MatrixType::random_uniform<RandomSource>(8);
        decomp(mat, 1E-11);
    }
    /* Complex case */ {
        using MatrixType = DenseMatrix<Complex<RealType>, MatrixMajor::Col, 3, 3>;
        const MatrixType temp{{{2, 1}, {-3, 6}, {12, 7}}, {{-50, -9}, {2, 180}, {-9, -6}}, {{-7, 8}, {546, 0}, {0, -25}}};
        const MatrixType mat = temp + temp.hermite();
        decomp(mat, 1E-12);
    }
    return 0;
}
