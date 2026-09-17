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

    void equality() {
        VectorND<float64> a = {1, 2, 3, 4};
        VectorND<float64> b = a;
        VectorND<float64> c = {1, 2, 3, 5};
        expect(a == b);
        expect(a != c);
        expect(a != b.head<2>());
        expect(a.segment<2>(1, 3) == b.segment<2>(1, 3));

        a.reserve(64);
        expect(a == b);

        VectorND<cfloat64> x = {cfloat64(1.0, 2.0), cfloat64(3.0, 4.0)};
        VectorND<cfloat64> y = x;
        VectorND<cfloat64> z = {cfloat64(1.0, 2.0), cfloat64(3.0, 5.0)};
        expect(x == y);
        expect(x != z);
    }
}

int main() {
    read();
    equality();
    return 0;
}
