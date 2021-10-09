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

#include "EqualityQuadraticProgramming.h"

namespace Physica::Core {
    /**
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006.448-496
     */
    template<class ScalarType>
    class QuadraticProgramming {
        DenseSymmMatrix<ScalarType, Dynamic> objectiveMatG;
        Vector<ScalarType, Dynamic> objectiveVecC;
        DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector> equalityConstraint;
        DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector> inequalityConstraint;
        Vector<ScalarType, Dynamic> x;
        Utils::Array<bool, Dynamic> activeConstraintFlags;
    public:
        template<class MatrixType, class VectorType>
        QuadraticProgramming(const LValueMatrix<MatrixType>& objectiveMatG_,
                             const LValueVector<VectorType>& objectiveVecC_,
                             const LValueMatrix<MatrixType>& equalityConstraint_,
                             const LValueMatrix<MatrixType>& inequalityConstraint_,
                             const LValueVector<VectorType>& initial);
        QuadraticProgramming(const QuadraticProgramming&) = delete;
        QuadraticProgramming(QuadraticProgramming&&) noexcept = delete;
        ~QuadraticProgramming() = default;
        /* Operators */
        QuadraticProgramming& operator=(const QuadraticProgramming&) = delete;
        QuadraticProgramming& operator=(QuadraticProgramming&&) noexcept = delete;
        /* Operations */
        void compute();
        /* Getters */
        [[nodiscard]] const Vector<ScalarType, Dynamic>& getSolution() const noexcept { return x; }
    };

    template<class ScalarType>
    template<class MatrixType, class VectorType>
    QuadraticProgramming<ScalarType>::QuadraticProgramming(const LValueMatrix<MatrixType>& objectiveMatG_,
                                                           const LValueVector<VectorType>& objectiveVecC_,
                                                           const LValueMatrix<MatrixType>& equalityConstraint_,
                                                           const LValueMatrix<MatrixType>& inequalityConstraint_,
                                                           const LValueVector<VectorType>& initial)
            : objectiveMatG(objectiveMatG_)
            , objectiveVecC(objectiveVecC_)
            , equalityConstraint(equalityConstraint_)
            , inequalityConstraint(inequalityConstraint_)
            , x(initial)
            , activeConstraintFlags(Vector<bool, Dynamic>(inequalityConstraint_.getRow(), false)) {
        assert(objectiveMatG.getRow() == objectiveVecC.getLength());
        assert(equalityConstraint.getColumn() == 0 || equalityConstraint.getColumn() == objectiveVecC.getLength() + 1);
        assert(inequalityConstraint.getColumn() == 0 || inequalityConstraint.getColumn() == objectiveVecC.getLength() + 1);
        compute();
    }

    template<class ScalarType>
    void QuadraticProgramming<ScalarType>::compute() {
        while (true) {

        };
    }
}
