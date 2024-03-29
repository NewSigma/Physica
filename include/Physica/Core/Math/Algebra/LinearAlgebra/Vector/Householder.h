/*
 * Copyright 2020-2022 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Conjugate.h"

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
    typename AnyVector::ScalarType::RealType householder(const RValueVector<AnyVector>& source,
                                                         LValueVector<OtherVector>& target) {
        using ScalarType = typename AnyVector::ScalarType;
        using RealType = typename ScalarType::RealType;
        assert(source.getLength() == target.getLength());

        const ScalarType v0 = source.calc(0);
        const RealType sourceNorm0 = v0.squaredNorm();
        const RealType squaredTailNorm = source.tail(1).squaredNorm();
        if (squaredTailNorm > std::numeric_limits<ScalarType>::min()) {
            const RealType norm = sqrt(squaredTailNorm + sourceNorm0);
            const ScalarType factor = v0.unit() * norm;
            const ScalarType factor1 = v0 + factor;
            const ScalarType factor2 = reciprocal(factor1);

            target.tail(1) = source.tail(1) * factor2;
            target[0] = (factor1 / factor).getReal();
            return norm;
        }
        target[0] = ScalarType::Two();
        target.tail(1) = ScalarType::Zero();
        return sqrt(sourceNorm0);
    }

    template<class AnyVector>
    typename AnyVector::ScalarType::RealType householderInPlace(LValueVector<AnyVector>& v) {
        return householder(v, v);
    }

    template<class MatrixType, class VectorType>
    void applyHouseholder(const RValueVector<VectorType>& householder, LValueMatrix<MatrixType>& mat) {
        using std::swap;
        using ScalarType = typename MatrixType::ScalarType;
        using T = Vector<ScalarType, VectorType::SizeAtCompile, VectorType::MaxSizeAtCompile>;
        T copy = householder;
        ScalarType temp = ScalarType::One();
        swap(temp, copy[0]);
        const T temp1 = copy * temp;
        mat -= temp1 * (copy.transpose().conjugate() * mat).compute();
    }

    template<class MatrixType, class VectorType>
    void applyHouseholder(LValueMatrix<MatrixType>& mat, const RValueVector<VectorType>& householder) {
        using std::swap;
        using ScalarType = typename MatrixType::ScalarType;
        using T = Vector<ScalarType, VectorType::SizeAtCompile, VectorType::MaxSizeAtCompile>;
        T copy = householder;
        ScalarType temp = ScalarType::One();
        swap(temp, copy[0]);
        using ProductType = decltype(mat * copy);
        using T1 = Vector<ScalarType, ProductType::SizeAtCompile, ProductType::MaxSizeAtCompile>;
        mat -= T1(mat * copy) * (copy.conjugate().transpose() * temp);
    }
}
