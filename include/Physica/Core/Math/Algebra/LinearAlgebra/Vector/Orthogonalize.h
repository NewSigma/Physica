/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/RValueMatrix.h"
#include "LValueVector.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.254-255
     */
    template<class MatrixType, class VectorType>
    void gramSchmidt(const RValueMatrix<MatrixType>& base, LValueVector<VectorType>& v) {
        assert(base.getRow() > base.getColumn());
        for (size_t i = 0; i < base.getColumn(); ++i) {
            auto col = base.col(i);
            auto dot = col.asVector().conjugate() * v.getDerived();
            v -= dot * col.asVector();
        }
        v.toUnit();
    }

    template<class MatrixType>
    void gramSchmidt(LValueMatrix<MatrixType>& m) {
        assert(m.getRow() >= m.getColumn());
        for (size_t i = 0; i < m.getColumn(); ++i) {
            auto col1 = m.col(i);
            col1.toUnit();
            for (size_t j = i + 1; j < m.getColumn(); ++j) {
                auto col2 = m.col(j);
                const auto dot = col1.asVector().conjugate() * col2.asVector();
                col2 -= dot * col1.asVector();
            }
        }
    }
}
