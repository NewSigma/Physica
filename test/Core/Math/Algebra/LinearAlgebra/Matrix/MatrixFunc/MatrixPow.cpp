/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixPow.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
template<Scalar T>
using Matrix2D = DenseMatrix<T, MatrixMajor::Col, 2, 2>;

namespace {
    template<Scalar T>
    Matrix2D<T> rotation(T theta) noexcept {
        Matrix2D<T> result{};
        result[0, 0] = cos(theta);
        result[0, 1] = sin(theta);
        result[1, 0] = -sin(theta);
        result[1, 1] = cos(theta);
        return result;
    }

    void gemv(int power) {
        using T = float64;
        const T theta = T::random_uniform<Random<PCG64DXSM, 1234>>();
        auto x = rotation(theta);
        Vector2D<T> proj{1, -1};
        Vector2D<T> answer = rotation(theta * T(power)) * proj;
        Vector2D<T> result = pow(x, power) * proj;
        expect(vectorNear(result, answer, 1E-15));
    }
}

int main() {
    gemv(2);
    gemv(3);
    return 0;
}
