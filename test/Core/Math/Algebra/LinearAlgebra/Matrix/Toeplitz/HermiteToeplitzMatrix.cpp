/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Toeplitz/HermiteToeplitzMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Test.h"

using namespace Physica;
using RandomSource = Random<>;

namespace {
    void assign() noexcept {
        using Tc = cfloat32;
        constexpr int N = 8;
        auto col = VectorND<Tc>::random_uniform<RandomSource>(N);
        col[0].imag() = 0;
        const auto toe = HermiteToeplitzMatrix<Tc>(std::move(col));
        MatrixND<Tc> buffer = toe;
        expect(buffer.isHermite());

        for (int i = 0; i < N - 1; ++i)
            for (int j = 1; j < N - i; ++j)
                expect(toe.getFirstCol()[i] == buffer[i + j, j]);
    }
}

int main() {
    static_assert(HermiteToeplitzMatrix<cfloat32>::isStaticHermite());
    assign();
    return 0;
}
