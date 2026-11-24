/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Optimization/SteepestDescent.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;
using Optimizer = SteepestDescent<ScalarType, 2>;
using VectorType = Optimizer::VectorType;

int main() {
    const auto func = [](VectorType x) -> ScalarType { return x.squaredNorm(); };
    const auto grad = [](VectorType x) -> VectorType { return x * ScalarType(2); };
    Optimizer sd(0.5, 1E-4, 0.1);
    sd.init({5, 4}, func, grad);
    do {
        sd.step(func, grad);
    } while(sd.getGradG().norm() > 0.001);
    expect(sd.getObjectiveValue() < 1E-8);
    return 0;
}
