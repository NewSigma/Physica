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
#include "Physica/Utils/Cycler.h"
#include "Physica/Utils/CUDA/DeviceProp.cuh"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Parallel/StreamPool.cuh"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;

void run(unsigned int sys, MatrixType& record, std::mt19937& gen);

int main() {
    Cycler::init();
    {
        MatrixType record(50, 1);
        std::mt19937::result_type seed;
        Physica::Utils::Random::rdrand(seed);
        std::mt19937 gen(seed);

        auto timeuse = Benchmark::run([&]() {
            run(0, record, gen);
        }, 8, 10);
        std::cout << "Time use(second): " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    return 0;
}
