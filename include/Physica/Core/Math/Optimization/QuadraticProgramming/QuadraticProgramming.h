/*
 * Copyright 2021-2022 WeiBo He.
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

#include <algorithm>
#include "EqualityQuadraticProgramming.h"

namespace Physica::Core {
    /**
     * Solve quadratic programming, that is
     * 
     * min 1/2 x^T G x + c^T x
     * s.t. A x = b
     *      C x >= d
     * 
     * equalityConstraint is matrix [A b]
     * inequalityConstraint is matrix [C d]
     * 
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006.448-496
     */
    template<class ScalarType>
    class QuadraticProgramming {
        using ConstraintMatrix = DenseMatrix<ScalarType, DenseMatrixOption::Row | DenseMatrixOption::Vector>;

        DenseSymmMatrix<ScalarType, Dynamic> objectiveMatG;
        Vector<ScalarType, Dynamic> objectiveVecC;
        ConstraintMatrix equalityConstraint;
        ConstraintMatrix inequalityConstraint;
        Vector<ScalarType, Dynamic> x;
        Utils::Array<bool, Dynamic> activeConstraintFlags;
    public:
        template<class MatrixType1, class VectorType1, class MatrixType2, class MatrixType3, class VectorType2>
        QuadraticProgramming(const LValueMatrix<MatrixType1>& objectiveMatG_,
                             const LValueVector<VectorType1>& objectiveVecC_,
                             const LValueMatrix<MatrixType2>& equalityConstraint_,
                             const LValueMatrix<MatrixType3>& inequalityConstraint_,
                             const LValueVector<VectorType2>& initial);
        QuadraticProgramming(const QuadraticProgramming&) = delete;
        QuadraticProgramming(QuadraticProgramming&&) noexcept = delete;
        ~QuadraticProgramming() = default;
        /* Operators */
        QuadraticProgramming& operator=(const QuadraticProgramming&) = delete;
        QuadraticProgramming& operator=(QuadraticProgramming&&) noexcept = delete;
        /* Operations */
        void compute(const ScalarType& precision = std::numeric_limits<ScalarType>::epsilon());
        void compute_nonconvex(const ScalarType& precision = std::numeric_limits<ScalarType>::epsilon());
        /* Getters */
        [[nodiscard]] const Vector<ScalarType, Dynamic>& getSolution() const noexcept { return x; }
        [[nodiscard]] ScalarType getValue() const;
    private:
        void updateVariables(const Vector<ScalarType, Dynamic>& direction);
        [[nodiscard]] ScalarType nextStepSize(const Vector<ScalarType, Dynamic>& direction, size_t& blockedAt);
        void updateActiveConstraints(ConstraintMatrix& activeConstraints);
    };

    template<class ScalarType>
    template<class MatrixType1, class VectorType1, class MatrixType2, class MatrixType3, class VectorType2>
    QuadraticProgramming<ScalarType>::QuadraticProgramming(const LValueMatrix<MatrixType1>& objectiveMatG_,
                                                           const LValueVector<VectorType1>& objectiveVecC_,
                                                           const LValueMatrix<MatrixType2>& equalityConstraint_,
                                                           const LValueMatrix<MatrixType3>& inequalityConstraint_,
                                                           const LValueVector<VectorType2>& initial)
            : objectiveMatG(objectiveMatG_)
            , objectiveVecC(objectiveVecC_)
            , equalityConstraint(equalityConstraint_)
            , inequalityConstraint(inequalityConstraint_)
            , x(initial)
            , activeConstraintFlags(Utils::Array<bool, Dynamic>(inequalityConstraint_.getRow(), false)) {
        assert(objectiveMatG.getRow() == objectiveVecC.getLength());
        assert(equalityConstraint.getColumn() == 0 || equalityConstraint.getColumn() == objectiveVecC.getLength() + 1);
        assert(inequalityConstraint.getColumn() == 0 || inequalityConstraint.getColumn() == objectiveVecC.getLength() + 1);
    }

    template<class ScalarType>
    void QuadraticProgramming<ScalarType>::compute(const ScalarType& precision) {
        ConstraintMatrix activeConstraints = equalityConstraint;
        while (true) {
            const EqualityQuadraticProgramming<ScalarType> EQP(objectiveMatG, objectiveVecC, activeConstraints, x);
            const Vector<ScalarType> vec_p = EQP.getSolution() - x;
            if (vec_p.norm() <= x.norm() * precision) {
                const auto& multipliers = EQP.getMultipliers();
                const auto minimum_ite = std::min_element(multipliers.cbegin(), multipliers.cend());
                if (!(*minimum_ite).isNegative())
                    break;
                
                auto minimum_index = std::distance(multipliers.cbegin(), minimum_ite);
                for (auto ite = activeConstraintFlags.begin(); ite != activeConstraintFlags.end(); ++ite) {
                    if (*ite == true) {
                        if (minimum_index == 0) {
                            *ite = false;
                            break;
                        }
                        --minimum_index;
                    }
                }
            }
            else
                updateVariables(vec_p);
            updateActiveConstraints(activeConstraints);
        };
    }
    /**
     * Return a local solution only
     */
    template<class ScalarType>
    void QuadraticProgramming<ScalarType>::compute_nonconvex(const ScalarType& precision) {
        ConstraintMatrix activeConstraints = equalityConstraint;
        while (true) {
            const EqualityQuadraticProgramming<ScalarType> EQP(objectiveMatG, objectiveVecC, activeConstraints, x);
            const Vector<ScalarType> vec_p = EQP.getSolution() - x;
            if (vec_p.norm() <= x.norm() * precision)
                break;
            else
                updateVariables(vec_p);
            updateActiveConstraints(activeConstraints);
        };
    }

    template<class ScalarType>
    ScalarType QuadraticProgramming<ScalarType>::getValue() const {
        return (x.copyToRowMatrix() * objectiveMatG).compute().row(0) * x * ScalarType(0.5) + objectiveVecC * x;
    }

    template<class ScalarType>
    void QuadraticProgramming<ScalarType>::updateVariables(const Vector<ScalarType, Dynamic>& direction) {
        size_t blockedAt;
        const ScalarType step = nextStepSize(direction, blockedAt);
        x = x + step * direction;
        if (step != ScalarType::One()) {
            assert(activeConstraintFlags[blockedAt] == false);
            activeConstraintFlags[blockedAt] = true;
        }
    }

    template<class ScalarType>
    ScalarType QuadraticProgramming<ScalarType>::nextStepSize(const Vector<ScalarType, Dynamic>& direction, size_t& blockedAt) {
        ScalarType result = ScalarType::One();
        for (size_t i = 0; i < activeConstraintFlags.getLength(); ++i) {
            const bool isActive = activeConstraintFlags[i];
            if (!isActive) {
                const auto row = inequalityConstraint.row(i);
                const auto head = row.head(row.getLength() - 1);
                const ScalarType dot = head * direction;
                if (dot.isNegative()) {
                    const ScalarType stepSize_i = (row[row.getLength() - 1] - head * x) / dot;
                    const bool less = stepSize_i < result;
                    if (less) {
                        result = stepSize_i;
                        blockedAt = i;
                    }
                }
            }
        }
        return result;
    }

    template<class ScalarType>
    void QuadraticProgramming<ScalarType>::updateActiveConstraints(ConstraintMatrix& activeConstraints) {
        const size_t activeInequality = std::count(activeConstraintFlags.begin(), activeConstraintFlags.end(), true);
        activeConstraints.resize(equalityConstraint.getRow() + activeInequality, x.getLength() + 1);

        size_t activeInequalityIndex = 0;
        for (size_t i = 0; i < activeConstraintFlags.getLength(); ++i) {
            if (activeConstraintFlags[i]) {
                auto row = activeConstraints.row(equalityConstraint.getRow() + activeInequalityIndex);
                row.asVector() = inequalityConstraint.row(i);
                ++activeInequalityIndex;
            }
        }
    }
}
