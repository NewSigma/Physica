/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void read() {
        using Tr = float32;
        using Tc = cfloat32;
        const auto complex = MatrixND<Tc>::random_uniform<RandomSource>(2, 3);

        VectorND<Tr> raw(complex.getSize() * 2, 0);
        raw.read(complex);

        MatrixND<Tc> result(2, 3, Tc(0));
        result.read(raw);
        expect(result == complex);
    }
}

int main() {
    read();
    return 0;
}
