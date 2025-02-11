/*
 * Copyright 2022-2024 Weibo He.
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
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/RValueMatrix.h"
#include "VectorImpl/LValueVector.h"

namespace Physica {
    /**
     * Reference:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:254-255
     */
    template<Matrix T, Vector U>
    void gramSchmidt(const T& base, U& v) {
        assert(base.getRow() > base.getCol() && "[Error]: base is over complete");
        for (size_t i = 0; i < base.getCol(); ++i) {
            const auto col = base.col(i);
            const auto dot = col.conjugate() * v;
            v -= dot * col;
        }
        v.toUnit();
    }
    /**
     * normGramSchmidt can run faster than gramSchmidt if \param base_ are normalized.
     * 
     * This method is vulnerable to numerical roundness.
     */
    template<Matrix T, Vector U>
    void normGramSchmidt(const T& base, U& v, typename U::ScalarType::RealType squaredNorm = 1) {
        using ScalarType = U::ScalarType;
        using RealType = ScalarType::RealType;
        [[maybe_unused]] constexpr double epsilon = 1E-3;
        assert(base.getRow() > base.getCol() && "[Error]: base is over complete");
        assert(scalarNear(v.squaredNorm(), squaredNorm, epsilon) && "[Error]: Invalid param");
        for (size_t i = 0; i < base.getCol(); ++i) {
            const auto col = base.col(i);
            [[maybe_unused]] const RealType temp = col.squaredNorm();
            assert(scalarNear(temp, RealType(1), epsilon) && "[Error]: Invalid param");
            const auto dot = col.conjugate() * v;
            v -= dot * col;
            squaredNorm -= dot.squaredNorm();
        }
        v *= reciprocal(sqrt(squaredNorm));
    }

    template<Matrix T>
    void gramSchmidt(T& m) {
        assert(m.getRow() >= m.getCol() && "[Error]: base is over complete");
        for (size_t i = 0; i < m.getCol(); ++i) {
            auto col1 = m.col(i);
            col1.toUnit();
            for (size_t j = i + 1; j < m.getCol(); ++j) {
                auto col2 = m.col(j);
                const auto dot = col1.conjugate() * col2;
                col2 -= dot * col1;
            }
        }
    }
}
