/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Statistics/LinearFit.h"
#include "Test.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using Fit = LinearFit<ScalarType>;

int main() {
    const auto x = VectorType::linspace(0, 4 * M_PI, 10);
    const VectorType y = sin(x);
    const VectorType p = Fit::polyfit(x, y, 7);
    const VectorType y1 = Fit::polyval(x, p);
    const ScalarType err = square(y - y1).sum() / ScalarType(y.getLength());
    expect(err < ScalarType(1E-3));
    return 0;
}
