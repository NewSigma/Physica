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
#pragma once

#include <cassert>
#include <iostream>
#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/MultiPrecision/ComplexScalar.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/LValueMatrix.h"

namespace Physica::Utils {
    template<class ScalarType>
    bool scalarNear(const Physica::Core::ScalarBase<ScalarType>& scalar1,
                    const Physica::Core::ScalarBase<ScalarType>& scalar2,
                    double precision) {
        assert(precision > 0);
        const auto& s1 = scalar1.getDerived();
        const auto& s2 = scalar2.getDerived();
        const ScalarType min = std::numeric_limits<ScalarType>::min();
        const bool useAbsCompare = (abs(s1) < min) || (abs(s2) < min);
        const ScalarType delta = s1 - s2;
        const ScalarType error = useAbsCompare ? abs(delta) : abs(delta / (s1 + s2) * ScalarType::Two());
        return error < ScalarType(precision);
    }

    template<class ScalarType>
    bool scalarNear(const Physica::Core::ComplexScalar<ScalarType>& s1,
                    const Physica::Core::ComplexScalar<ScalarType>& s2,
                    double precision) {
        assert(precision > 0);
        return scalarNear((s1 - s2).norm(), typename ScalarType::RealType(0), precision);
    }

    template<class VectorType1, class VectorType2>
    bool vectorNear(const Physica::Core::LValueVector<VectorType1>& v1,
                    const Physica::Core::LValueVector<VectorType2>& v2,
                    double precision) {
        assert(v1.getLength() == v2.getLength());
        for (size_t i = 0; i < v1.getLength(); ++i)
            if (!scalarNear(v1[i], v2[i], precision))
                return false;
        return true;
    }

    template<class MatrixType, class MatrixType2>
    bool matrixNear(const Physica::Core::LValueMatrix<MatrixType>& m1,
                    const Physica::Core::LValueMatrix<MatrixType2>& m2,
                    double precision) {
        assert(m1.getRow() == m2.getRow());
        assert(m1.getColumn() == m2.getColumn());
        for (size_t i = 0; i < m1.getColumn(); ++i)
            for (size_t j = 0; j < m1.getRow(); ++j)
                if (!scalarNear(m1(j, i), m2(j, i), precision))
                    return false;
        return true;
    }
}
