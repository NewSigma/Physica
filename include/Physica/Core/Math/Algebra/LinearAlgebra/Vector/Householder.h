/*
 * Copyright 2020-2021 WeiBo He.
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

#include "LValueVector.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/LValueMatrix.h"

namespace Physica::Core {
    /**
     * The first element of \param target will be the factor to construct houseHolder matrix.
     * The other parts of \param target will be essential HouseHolder vector.
     * 
     * \return The norm of \param source
     * 
     * References:
     * [1] Golub, GeneH. Matrix computations = 矩阵计算 / 4th edition[M]. 人民邮电出版社, 2014.
     * [2] Eigen https://eigen.tuxfamily.org/
     */
    template<class AnyVector, class OtherVector>
    typename AnyVector::ScalarType householder(const LValueVector<AnyVector>& __restrict source,
                                               LValueVector<OtherVector>& __restrict target) {
        using ScalarType = typename AnyVector::ScalarType;
        assert(source.getLength() == target.getLength());
        const ScalarType abs_first = abs(source[0]);
        const bool minus = source[0].isNegative();
        const ScalarType norm = source.getDerived().norm();
        ScalarType factor = reciprocal(ScalarType(abs_first + norm));
        if (minus)
            factor.toOpposite();

        target.tail(1) = source.tail(1) * factor;
        target[0] = ScalarType(1) + abs_first / norm;
        return norm;
    }

    template<class AnyVector>
    typename AnyVector::ScalarType householderInPlace(LValueVector<AnyVector>& v) {
        return householder(v, v);
    }

    template<class MatrixType, class VectorType>
    void applyHouseholder(const LValueVector<VectorType>& householder, LValueMatrix<MatrixType>& mat) {
        using ScalarType = typename MatrixType::ScalarType;
        Vector<ScalarType> copy = householder;
        ScalarType temp = ScalarType::One();
        std::swap(temp, copy[0]);
        const auto mat1 = (copy * temp).copyToColMatrix();
        mat -= mat1 * (copy.copyToRowMatrix() * mat).compute();
    }

    template<class MatrixType, class VectorType>
    void applyHouseholder(LValueMatrix<MatrixType>& mat, const LValueVector<VectorType>& householder) {
        using ScalarType = typename MatrixType::ScalarType;
        Vector<ScalarType> copy = householder;
        ScalarType temp = ScalarType::One();
        std::swap(temp, copy[0]);
        const auto mat1 = (copy * temp).copyToRowMatrix();
        mat -= (mat * copy.copyToColMatrix()).compute() * mat1;
    }
}
