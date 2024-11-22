/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica::Core;
using ScalarType = float64;
using VectorType = Diff<VectorND<ScalarType>, DiffMode::Reverse, 1>;

int main() {
    VectorType v = VectorType::random_uniform(16, Random<MT19937, std::mt19937::default_seed>::getInstance());
    auto sum = v.sum();
    sum.reverse();
    for (size_t i = 0; i < v.getLength(); ++i)
        if (v.calc(i).getGrad() != ScalarType(1))
            return 1;
    return 0;
}
