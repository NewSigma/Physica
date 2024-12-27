/*
 * Copyright 2020-2024 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/LValueMatrix.h"

namespace Physica::Core {
    /**
     * The first element of \param target will be the factor to construct houseHolder matrix.
     * The other parts of \param target will be essential HouseHolder vector.
     * 
     * \return The norm of \param source
     * 
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013
     * [2] Eigen; https://eigen.tuxfamily.org/
     */
    template<Vector T, LVector U>
    T::ScalarType::RealType householder(const T& source, U& target) {
        using ScalarType = T::ScalarType;
        using RealType = ScalarType::RealType;
        assert(source.getLength() == target.getLength());

        const ScalarType v0 = source.calc(0);
        const RealType sourceNorm0 = v0.squaredNorm();
        const RealType squaredTailNorm = source.tail(1).squaredNorm();
        if (squaredTailNorm > std::numeric_limits<ScalarType>::min()) [[likely]] {
            const RealType norm = sqrt(squaredTailNorm + sourceNorm0);
            const ScalarType factor = v0.unit() * norm;
            const ScalarType factor1 = v0 + factor;
            const ScalarType factor2 = reciprocal(factor1);

            target.tail(1) = source.tail(1) * factor2;
            target[0] = (factor1 / factor).real();
            return norm;
        }
        else {
            const bool isZeroVector = sourceNorm0 < std::numeric_limits<ScalarType>::min();
            if (isZeroVector) {
                target = RealType(0);
                return RealType(0);
            }
        }
        target[0] = RealType(2);
        target.tail(1) = RealType(0);
        return sqrt(sourceNorm0);
    }

    template<LVector T>
    T::ScalarType::RealType householderInPlace(T& v) {
        return householder(v, v);
    }

    template<LMatrix M, Vector V>
    void applyHouseholder(const V& householder, M& mat) {
        using ScalarType = M::ScalarType;
        using BufferType = DenseVector<ScalarType, V::SizeAtCompile>;
        BufferType copy = householder;
        auto temp = ScalarType(1);
        temp.swap(copy[0]);
        const BufferType temp1 = copy * temp;
        mat -= temp1 * (copy.hermite() * mat).compute();
    }

    template<LMatrix M, Vector V>
    void applyHouseholder(M& mat, const V& householder) {
        using ScalarType = M::ScalarType;
        using BufferType = DenseVector<ScalarType, V::SizeAtCompile>;
        BufferType copy = householder;
        ScalarType temp = ScalarType(1);
        temp.swap(copy[0]);
        using ProductType = decltype(mat * copy);
        using BufferType1 = DenseVector<ScalarType, ProductType::SizeAtCompile>;
        mat -= BufferType1(mat * copy) * (copy.hermite() * temp);
    }
}
