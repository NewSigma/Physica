/*
 * Copyright 2021 WeiBo He.
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
#include <cassert>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/LValueMatrix.h"

template<class ScalarType>
bool floatNear(const ScalarType& s1, const ScalarType& s2, double precision) {
    assert(precision > 0);
    return s2.isZero() ? (abs(s1) < ScalarType(precision)) : (abs((s1 - s2) / s2) < ScalarType(precision));
}

template<class MatrixType, class MatrixType2>
bool matrixNear(const Physica::Core::LValueMatrix<MatrixType>& m1,
                const Physica::Core::LValueMatrix<MatrixType2>& m2,
                double precision) {
    assert(m1.getRow() == m2.getRow());
    assert(m1.getColumn() == m2.getColumn());
    for (size_t i = 0; i < m1.getColumn(); ++i)
        for (size_t j = 0; j < m1.getRow(); ++j)
            if (!floatNear(m1(j, i), m2(j, i), precision))
                return false;
    return true;
}
