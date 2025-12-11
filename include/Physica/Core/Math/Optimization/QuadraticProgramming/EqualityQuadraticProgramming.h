/*
 * Copyright 2021-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/DenseLU.h"

namespace Physica {
    /**
     * Solve quadratic programming with equality constraints only, that is
     * 
     * min 1/2 x^T G x + c^T x
     * s.t. A x = b
     * 
     * constraints is matrix [A b]
     * 
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006:448-496
     */
    template<Scalar T>
    class EqualityQuadraticProgramming {
        DenseSymmMatrix<T> objectiveMatG;
        VectorND<T> objectiveVecC;
        DenseMatrix<T, MatrixOption::Row> constraints;
        VectorND<T> x;
        VectorND<T> multipliers;
    public:
        EqualityQuadraticProgramming(const Matrix auto& objectiveMatG_,
                                     const Vector auto& objectiveVecC_,
                                     const Matrix auto& constraints_,
                                     const Vector auto& initial);
        EqualityQuadraticProgramming(const EqualityQuadraticProgramming&) = delete;
        EqualityQuadraticProgramming(EqualityQuadraticProgramming&&) noexcept = delete;
        ~EqualityQuadraticProgramming() = default;
        /* Operators */
        EqualityQuadraticProgramming& operator=(const EqualityQuadraticProgramming&) = delete;
        EqualityQuadraticProgramming& operator=(EqualityQuadraticProgramming&&) noexcept = delete;
        /* Operations */
        void compute();
        /* Getters */
        [[nodiscard]] const VectorND<T>& getSolution() const noexcept { return x; }
        [[nodiscard]] const VectorND<T>& getMultipliers() const noexcept { return multipliers; }
    };

    template<Scalar T>
    EqualityQuadraticProgramming<T>::EqualityQuadraticProgramming(const Matrix auto& objectiveMatG_,
                                                                  const Vector auto& objectiveVecC_,
                                                                  const Matrix auto& constraints_,
                                                                  const Vector auto& initial)
            : objectiveMatG(objectiveMatG_)
            , objectiveVecC(objectiveVecC_)
            , x(initial) {
        assert(objectiveMatG.getRow() == objectiveVecC.getLength());
        assert(constraints_.getCol() == 0 || constraints_.getCol() == objectiveVecC.getLength() + 1);
        assert(x.getLength() == objectiveVecC.getLength());
        if (constraints_.getSize() != 0)
            constraints = constraints_;
        compute();
    }

    template<Scalar T>
    void EqualityQuadraticProgramming<T>::compute() {
        const size_t degreeOfFreedom = x.getLength();
        const bool haveConstraints = constraints.getRow() != 0;
        VectorND<T> equationVecB(objectiveMatG.getRow() + constraints.getRow());
        /* Assemble vector */ {
            const VectorND<T> g = objectiveMatG * x + objectiveVecC;
            auto head = equationVecB.head(degreeOfFreedom);
            head = g;

            if (haveConstraints) {
                const auto matA = constraints.leftCols(degreeOfFreedom);
                const auto vecB = constraints.col(degreeOfFreedom);
                const VectorND<T> h = matA * x - vecB;
                auto tail = equationVecB.tail(degreeOfFreedom);
                tail = h;
            }
        }

        const size_t problemSize = degreeOfFreedom + constraints.getRow();
        DenseMatrix<T, MatrixOption::Row> equationMatA(problemSize, problemSize);
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
                bottomRight = T(0);
            }
            else
                equationMatA = objectiveMatG;
        }
        const DenseMatrix<T, MatrixOption::Row> inv_equationMatA = equationMatA.inv();
        const VectorND<T> solution = inv_equationMatA * equationVecB;
        x -= solution.head(degreeOfFreedom);
        if (haveConstraints)
            multipliers = solution.tail(degreeOfFreedom);
    }
}