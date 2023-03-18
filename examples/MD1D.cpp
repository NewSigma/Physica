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
#include "Physica/Core/Physics/MD/MD1D.h"
#include "Physica/Utils/Random.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double, false>;
constexpr double timeStep = 0.001;

int main() {
    std::mt19937::result_type seed;
    Physica::Utils::Random::rdrand(seed);
    std::mt19937 gen(seed);

    MD1D<ScalarType> md(200, 1024, timeStep, timeStep * 0.01);
    md.init(gen);
    for (int i = 0; i < 100000; ++i)
        md.nve_step();
    return 0;
}
