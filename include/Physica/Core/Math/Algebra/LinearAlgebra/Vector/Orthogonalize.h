/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
