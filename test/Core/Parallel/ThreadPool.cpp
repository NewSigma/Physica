/*
 * Copyright 2021-2025 Weibo He.
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
#include "Physica/Core/Parallel/Executor/SeqExecutor.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"

using namespace Physica;
using ScalarType = float64;
using MatrixType = DenseMatrix<ScalarType>;
using VectorType = VectorND<ScalarType>;
using RandomSource = Random<MT19937>;

void func(size_t) {
    printf("Thread ID: %d\n", ThreadPool::getThreadID());
}

int main() {
    const MatrixType A = MatrixType::random_uniform<RandomSource>(4);
    ThreadPool::numThreadRequired = 4;
    ThreadPool& pool = ThreadPool::getInstance();

    VectorType result_seq(4);
    VectorType result_par(4);
    auto sum_col = [&](VectorType& result, unsigned int i) { result[i] = A.col(i).sum(); };
    SeqExecutor::parallel_for([=, &result_seq](unsigned int i) { sum_col(result_seq, i); }, 4, 4);
    ThreadExecutor::parallel_for([=, &result_par](unsigned int i) { sum_col(result_par, i); }, 4, 4).wait();

    pool.shouldExit();
    return result_seq != result_par;
}
