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
     * Solve quadratic programming with equality constraints only, that is
     * 
     * min 1/2 x^T G x + c^T x
     * s.t. A x = b
     * 
     * constraints is matrix [A b]
     * 
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
        assert(constraints_.getColumn() == 0 || constraints_.getColumn() == objectiveVecC.getLength() + 1);
        assert(x.getLength() == objectiveVecC.getLength());
        compute();
    }

    template<class ScalarType>
    void EqualityQuadraticProgramming<ScalarType>::compute() {
        const size_t degreeOfFreedom = x.getLength();
        const bool haveConstraints = constraints.getRow() != 0;
        Vector<ScalarType, Dynamic> equationVecB(objectiveMatG.getRow() + constraints.getRow());
        /* Assemble vector */ {
            const Vector<ScalarType, Dynamic> g = (objectiveMatG * x.copyToColMatrix()).compute().col(0) + objectiveVecC;
            auto head = equationVecB.head(degreeOfFreedom);
            head = g;

            if (haveConstraints) {
                const auto matA = constraints.leftCols(degreeOfFreedom);
                const auto vecB = constraints.col(degreeOfFreedom);
                const Vector<ScalarType, Dynamic> h = (matA * x.copyToColMatrix()).compute().col(0) - vecB.asVector();
                auto tail = equationVecB.tail(degreeOfFreedom);
                tail = h;
            }
        }

        const size_t problemSize = degreeOfFreedom + constraints.getRow();
        DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector> equationMatA(problemSize, problemSize);
        /* Assemble matrix */ {
            if (haveConstraints) {
                const auto matA = constraints.leftCols(degreeOfFreedom);
                auto topLeft = equationMatA.topLeftCorner(degreeOfFreedom);
                topLeft = objectiveMatG;
                auto bottomLeft = equationMatA.bottomLeftCorner(degreeOfFreedom, degreeOfFreedom);
                bottomLeft = matA;
                auto topRight = equationMatA.topRightCorner(degreeOfFreedom, degreeOfFreedom);
                topRight = matA.transpose();
                auto bottomRight = equationMatA.bottomRightCorner(degreeOfFreedom);
                bottomRight = ScalarType::Zero();
            }
            else
                equationMatA = objectiveMatG;
        }
        const DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector> inv_equationMatA = equationMatA.inverse();
        const auto solution = (inv_equationMatA * equationVecB.moveToColMatrix()).compute();
        const auto col = solution.col(0);
        x -= col.head(degreeOfFreedom);
        if (haveConstraints)
            multipliers = col.tail(degreeOfFreedom);
    }
}