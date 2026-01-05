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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixFunction/MatrixPow.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
template<Scalar T>
using Matrix2D = DenseMatrix<T, MatrixOption::Col, 2, 2>;

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
        const T theta = T::random_uniform<Random<>>();
        auto x = rotation(theta);
        Vector2D<T> proj{1, -1};
        Vector2D<T> answer = rotation(theta * T(power)) * proj;
        Vector2D<T> result = pow(x, power) * proj;
        expect(vectorNear(result, answer, 1E-14));
    }
}

int main() {
    gemv(2);
    gemv(3);
    return 0;
}
