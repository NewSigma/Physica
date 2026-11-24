/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DiffDenseMatrix.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<>;

namespace {
    void sum() {
        using M = DenseMatrix<Diff<T, DiffMode::Reverse, 1>>;
        M m = M::random_uniform<RandomSource>(4, 4);
        m.sum().reverse();

        auto v = m.flatten();
        for (size_t i = 0; i < v.getLength(); ++i)
            expect(v.calc(i).grad() == T(1));
    }
}

int main() {
    sum();
    return 0;
}
