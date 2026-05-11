/*
 * Copyright 2022-2026 Weibo He.
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
     * Use snake case to respect the authors' contributions
     *
     * Reference:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:254-255
     */
    void gram_schmidt(const Matrix auto& base, Vector auto& v) {
        assert(base.getRow() > base.getCol() && "[Error]: base is over complete");
        for (size_t i = 0; i < base.getCol(); ++i) {
            const auto col = base.col(i);
            const auto dot = col.conjugate() * v;
            v -= dot * col;
        }
        v.toUnit();
    }

    void gram_schmidt(const Matrix auto& base, Vector auto& v, Vector auto&& dots) {
        assert(base.getRow() > base.getCol() && "[Error]: base is over complete");
        assert(dots.getLength() == base.getCol());
        for (size_t i = 0; i < base.getCol(); ++i) {
            const auto col = base.col(i);
            const auto dot = col.conjugate() * v;
            v -= dot * col;
            dots[i] = dot;
        }
        v.toUnit();
    }

    void gram_schmidt(Matrix auto& m) {
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
