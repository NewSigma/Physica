/*
 * Copyright 2021 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixChain.h"

using namespace Physica;
using T = float64;
using MatrixType = DenseMatrix<T>;
using RandomSource = Random<MT19937, 10000>;

int main() {
    MatrixType a[6];
    a[0] = MatrixType::random_uniform<RandomSource>(30, 35);
    a[1] = MatrixType::random_uniform<RandomSource>(35, 15);
    a[2] = MatrixType::random_uniform<RandomSource>(15, 5);
    a[3] = MatrixType::random_uniform<RandomSource>(5, 10);
    a[4] = MatrixType::random_uniform<RandomSource>(10, 20);
    a[5] = MatrixType::random_uniform<RandomSource>(20, 25);
    /**
    * Reference:
    * [1] Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, Clifford Stein . 算法导论(第三版)[M]. 北京: 机械工业出版社, 2013:215
    */
    const MatrixType answer = (a[0] * (a[1] * a[2]).compute()).compute() * ((a[3] * a[4]).compute() * a[5]).compute();

    auto chain = MatrixChain<T>(6);
    for (int i = 0; i < 6; ++i)
        chain[i] = std::move(a[i]);
    chain.dynamicProgram();
    const auto result = chain.multiply();
    if (!matrixNear(result, answer, 1E-15))
        return 1;
    return 0;
}
