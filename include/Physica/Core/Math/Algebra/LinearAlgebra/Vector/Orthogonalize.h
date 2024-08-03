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

namespace Physica::Core {
    /**
     * Reference:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:254-255
     */
    template<class MatrixType, class VectorType>
    void gramSchmidt(const RValueMatrix<MatrixType>& base_, LValueVector<VectorType>& v) {
        const auto& base = base_.getDerived();
        assert(base.getRow() > base.getColumn() && "[Error]: base is over complete");
        for (size_t i = 0; i < base.getColumn(); ++i) {
            const auto col = base.col(i);
            const auto dot = col.asVector().conjugate() * v.getDerived();
            v -= dot * col.asVector();
        }
        v.toUnit();
    }
    /**
     * normGramSchmidt can run faster than gramSchmidt if \param base_ are normalized.
     * 
     * This method is vulnerable to numerical roundness.
     */
    template<class MatrixType, class VectorType>
    void normGramSchmidt(
            const RValueMatrix<MatrixType>& base_,
            LValueVector<VectorType>& v,
            typename VectorType::ScalarType::RealType squaredNorm = 1) {
        using ScalarType = typename VectorType::ScalarType;
        using RealType = typename ScalarType::RealType;
        [[maybe_unused]] constexpr double epsilon = 1E-3;
        const auto& base = base_.getDerived();
        assert(base.getRow() > base.getColumn() && "[Error]: base is over complete");
        assert(scalarNear(v.squaredNorm(), squaredNorm, epsilon) && "[Error]: Invalid param");
        for (size_t i = 0; i < base.getColumn(); ++i) {
            const auto col = base.col(i);
            [[maybe_unused]] const RealType temp = col.squaredNorm();
            assert(scalarNear(temp, RealType(1), epsilon) && "[Error]: Invalid param");
            const auto dot = col.asVector().conjugate() * v.getDerived();
            v -= dot * col.asVector();
            squaredNorm -= dot.squaredNorm();
        }
        v *= reciprocal(sqrt(squaredNorm));
    }

    template<class MatrixType>
    void gramSchmidt(LValueMatrix<MatrixType>& m_) {
        auto& m = m_.getDerived();
        assert(m.getRow() >= m.getColumn() && "[Error]: base is over complete");
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
