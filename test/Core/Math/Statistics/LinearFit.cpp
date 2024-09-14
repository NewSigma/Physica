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
#include <iostream>
#include <Physica/Core/Math/Statistics/LinearFit.h>

using namespace Physica::Core;
using ScalarType = float64;
using VectorType = Vector<ScalarType>;
using Fit = LinearFit<ScalarType>;

int main() {
    const auto x = VectorType::linspace(0, 4 * M_PI, 10);
    const VectorType y = sin(x);
    const VectorType p = Fit::polyfit(x, y, 7);
    const VectorType y1 = Fit::polyval(x, p);
    const ScalarType err = square(y - y1).sum() / ScalarType(y.getLength());
    if (err > ScalarType(1E-3))
        return 1;
    return 0;
}
