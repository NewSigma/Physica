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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Vector.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006.448-496
     */
    template<class ScalarType>
    class EqualityQuadraticProgramming {
        DenseSymmMatrix<ScalarType, Dynamic> objectiveMatG;
        Vector<ScalarType, Dynamic> objectiveVecC;
        DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector> constraints;
        Vector<ScalarType, Dynamic> x;
        Vector<ScalarType, Dynamic> multipliers;
    public:
        template<class MatrixType1, class VectorType1, class MatrixType2, class VectorType2>
        EqualityQuadraticProgramming(const RValueMatrix<MatrixType1>& objectiveMatG_,
                                     const RValueVector<VectorType1>& objectiveVecC_,
                                     const RValueMatrix<MatrixType2>& constraints_,
                                     const RValueVector<VectorType2>& initial);
        EqualityQuadraticProgramming(const EqualityQuadraticProgramming&) = delete;
        EqualityQuadraticProgramming(EqualityQuadraticProgramming&&) noexcept = delete;
        ~EqualityQuadraticProgramming() = default;
        /* Operators */
        EqualityQuadraticProgramming& operator=(const EqualityQuadraticProgramming&) = delete;
        EqualityQuadraticProgramming& operator=(EqualityQuadraticProgramming&&) noexcept = delete;
        /* Operations */
        void compute();
        /* Getters */
        [[nodiscard]] const Vector<ScalarType, Dynamic>& getSolution() const noexcept { return x; }
        [[nodiscard]] const Vector<ScalarType, Dynamic>& getMultipliers() const noexcept { return multipliers; }
    };

    template<class ScalarType>
    template<class MatrixType1, class VectorType1, class MatrixType2, class VectorType2>
    EqualityQuadraticProgramming<ScalarType>::EqualityQuadraticProgramming(const RValueMatrix<MatrixType1>& objectiveMatG_,
                                                                           const RValueVector<VectorType1>& objectiveVecC_,
                                                                           const RValueMatrix<MatrixType2>& constraints_,
                                                                           const RValueVector<VectorType2>& initial)
            : objectiveMatG(objectiveMatG_)
            , objectiveVecC(objectiveVecC_)
            , constraints(constraints_)
            , x(initial) {
        assert(objectiveMatG.getRow() == objectiveVecC.getLength());
        assert(constraints_.getColumn() == objectiveVecC.getLength() + 1);
        assert(x.getLength() == objectiveVecC.getLength());
        compute();
    }

    template<class ScalarType>
    void EqualityQuadraticProgramming<ScalarType>::compute() {
        const Vector<ScalarType, Dynamic> g = (objectiveMatG * x.copyToColMatrix()).compute().col(0) + objectiveVecC;
        
        const size_t length = x.getLength();
        const auto matA = constraints.leftCols(length);
        const auto vecB = constraints.col(length);
        const Vector<ScalarType, Dynamic> h = (matA * x.copyToColMatrix()).compute().col(0) - vecB.asVector();
        
        const size_t totalLength = g.getLength() + h.getLength();
        Vector<ScalarType, Dynamic> equationVecB(totalLength);
        /* Assemble vector */ {
            auto head = equationVecB.head(g.getLength());
            head = g;
            auto tail = equationVecB.tail(g.getLength());
            tail = h;
        }

        DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector> equationMatA(totalLength, totalLength);
        /* Assemble matrix */ {
            auto topLeft = equationMatA.topLeftCorner(g.getLength());
            topLeft = objectiveMatG;
            auto bottomLeft = equationMatA.bottomLeftCorner(g.getLength(), g.getLength());
            bottomLeft = matA;
            auto topRight = equationMatA.topRightCorner(g.getLength(), g.getLength());
            topRight = matA.transpose();
            auto bottomRight = equationMatA.bottomRightCorner(g.getLength());
            bottomRight = ScalarType::Zero();
        }
        const DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector> inv_equationMatA = equationMatA.inverse();
        const auto solution = (inv_equationMatA * equationVecB.moveToColMatrix()).compute();
        const auto col = solution.col(0);
        x -= col.head(length);
        multipliers = col.tail(length);
    }
}