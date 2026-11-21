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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/MatrixProduct/BatchedGEMM.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

int main() {
    using T = float32;
    const size_t batch = 3, row = 7, col = 5, inner = 6;
    const T alpha(1), beta(0);
    auto a = MatrixND<T>::template random_uniform<RandomSource>(row, batch * inner);
    auto b = MatrixND<T>::template random_uniform<RandomSource>(inner, batch * col);
    auto c = MatrixND<T>::template random_uniform<RandomSource>(row, batch * col);

    MatrixND<T> answer(row, batch * col);
    for (size_t s = 0; s < batch; ++s) {
        auto answerBlock = answer.block(0, row, s * col, col);
        answerBlock = a.block(0, row, s * inner, inner) * b.block(0, inner, s * col, col);
    }

    auto aFirst = a.block(0, row, 0, inner);
    auto bFirst = b.block(0, inner, 0, col);
    auto cFirst = c.block(0, row, 0, col);
    BatchedGEMM<T> gemm(aFirst, bFirst, cFirst, batch);
    for (size_t s = 0; s < batch; ++s)
        gemm.batch(s, {a.data() + s * row * inner, b.data() + s * inner * col, c.data() + s * row * col});

    if constexpr (HasMKL()) {
        gemm.compute(alpha, beta);
        expect<RandomSource>(matrixNear(c, answer, 1E-4));
    }
    return 0;
}
