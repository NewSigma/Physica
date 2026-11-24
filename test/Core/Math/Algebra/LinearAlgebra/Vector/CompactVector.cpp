/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
