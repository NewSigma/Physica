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
#include <iostream>
#include <gperftools/profiler.h>
#include "Physica/Utils/Random.h"
#include "Physica/Utils/BenchmarkHelper.h"

using namespace Physica::Utils;

void run(std::mt19937& gen);

int main() {
    Cycler::init();

    std::mt19937::result_type seed;
    Physica::Utils::Random::rdrand(seed);
    std::mt19937 gen(seed);

    auto pair = Benchmark::run([&gen]() {
        run(gen);
    }, 6, 6);
    std::cout << pair.first << '(' << pair.second << ')' << std::endl;
    return 0;
}
