/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.cuh"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using RandomSource = Random<>;
using T = float32;

int main() {
    using VectorType = device_obj<VectorND<Diff<T, DiffMode::Reverse>>>;
    auto v = VectorType(8, 0);
    v.reverse(T(1));
    CUDAContext::getInstance().wait();
    const auto grads = v.grads().toHost();

    for (auto elem : grads)
        if (elem != T(1))
            return 1;
    return 0;
}
