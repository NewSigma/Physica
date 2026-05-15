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
