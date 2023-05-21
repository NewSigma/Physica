/*
 * Copyright 2023 WeiBo He.
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
#include <random>
#include <fstream>
#include <iostream>
#include <gperftools/profiler.h>
#include "Physica/Utils/Random.h"
#include "Physica/Utils/BenchmarkHelper.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Parallel/StreamPool.cuh"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;

void run(unsigned int sys, MatrixType& record, std::mt19937& gen);

int main() {
    MatrixType record(1000, 8);
    VectorType mean(record.getRow()), devia(record.getRow());
    ThreadPool::numThreadRequired = 8;
    {
        ThreadExecutor::parallel_for([&record](unsigned int sys) {
            std::mt19937::result_type seed;
            Physica::Utils::Random::rdrand(seed);
            std::mt19937 gen(seed);

            run(sys, record, gen);
        }, record.getColumn(), ThreadPool::numThreadRequired).wait();
        ThreadPool::getInstance().shouldExit();

        for (size_t i = 0; i < mean.getLength(); ++i) {
            mean[i] = Physica::Core::mean(record.row(i));
            devia[i] = Physica::Core::deviation(record.row(i));
        }
        const ScalarType factor = reciprocal(mean[0]);
        mean *= factor;
        devia *= factor;

        std::ofstream fout("data");
        fout << mean << devia;
    }
    return 0;
}
