/*
 * Copyright 2023 Weibo He.
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
#include "Physica/Core/AI/NeuralNetwork/Loss.h"

using namespace Physica;
using ScalarType = float32;

void testSelectProperty() {
    const VectorND<ScalarType> result{-3.34036088, -109.5531235, 13.51656151, 11.29175949};
    const ScalarType l1 = Loss<ScalarType>::crossEntropy(result, 3);
    if (!std::isfinite(double(l1)))
        exit(EXIT_FAILURE);

    const ScalarType l2 = Loss<ScalarType>::crossEntropy(result, 1);
    if (!std::isfinite(double(l2)))
        exit(EXIT_FAILURE);
}

void testOverFlow() {
    const VectorND<ScalarType> result{555.321167, 364.9577942, 355.3863831, -594.8062134};
    for (size_t i = 0; i < result.getLength(); ++i) {
        const ScalarType s = Loss<ScalarType>::softmax(result, i);
        if (!std::isfinite(double(s)))
            exit(EXIT_FAILURE);
    }
}

int main() {
    testSelectProperty();
    testOverFlow();
    return 0;
}
