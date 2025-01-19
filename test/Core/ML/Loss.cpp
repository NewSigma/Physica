/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/ML/NeuralNetwork/Loss.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DiffVector.h"

using namespace Physica;
using RandomType = Random<MT19937, 10000>;

void testSelectProperty() {
    const VectorND<float32> result{-3.34036088, -109.5531235, 13.51656151, 11.29175949};
    const float32 l1 = crossEntropy(result, 3);
    if (!l1.isFinite())
        exit(EXIT_FAILURE);

    const float32 l2 = crossEntropy(result, 1);
    if (!l2.isFinite())
        exit(EXIT_FAILURE);
}

void testOverFlow() {
    const VectorND<float32> result{555.321167, 364.9577942, 355.3863831, -594.8062134};
    for (size_t i = 0; i < result.getLength(); ++i) {
        const float32 s = softmax(result, i);
        if (!s.isFinite())
            exit(EXIT_FAILURE);
    }
}

void testDiff() {
    using dfloat = Diff<float32, DiffMode::Reverse, 1>;
    constexpr int Label = 0;
    const auto x = VectorND<dfloat>::random_uniform<RandomType>(8);
    auto loss = [&x](size_t k) -> float32 {
        float32 result = 0;
        for (size_t i = 0; i < x.getLength(); ++i)
            result += softmax(x.values(), i) * (float32(i == k) - float32(Label == k));
        return result;
    };

    crossEntropy(x, Label).reverse();
    for (size_t i = 0; i < x.getLength(); ++i) {
        if (!scalarNear(x.grads()[i], loss(i), 1E-6))
            exit(EXIT_FAILURE);
    }
}

int main() {
    testSelectProperty();
    testOverFlow();
    testDiff();
    return 0;
}
