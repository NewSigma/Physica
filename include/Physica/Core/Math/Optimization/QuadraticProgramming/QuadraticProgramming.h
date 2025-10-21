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

#include <algorithm>
#include "EqualityQuadraticProgramming.h"

namespace Physica {
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
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006:448-496
     */
    template<Scalar T>
    class QuadraticProgramming {
        using ConstraintMatrix = DenseMatrix<T, MatrixOption::Row>;

        DenseSymmMatrix<T> objectiveMatG;
        VectorND<T> objectiveVecC;
        ConstraintMatrix equalityConstraint;
        ConstraintMatrix inequalityConstraint;
        VectorND<T> x;
        Array<bool> activeConstraintFlags;
    public:
        QuadraticProgramming(const Matrix auto& objectiveMatG_,
                             const Vector auto& objectiveVecC_,
                             const Matrix auto& equalityConstraint_,
                             const Matrix auto& inequalityConstraint_,
                             const Vector auto& initial);
        QuadraticProgramming(const QuadraticProgramming&) = delete;
        QuadraticProgramming(QuadraticProgramming&&) noexcept = delete;
        ~QuadraticProgramming() = default;
        /* Operators */
        QuadraticProgramming& operator=(const QuadraticProgramming&) = delete;
        QuadraticProgramming& operator=(QuadraticProgramming&&) noexcept = delete;
        /* Operations */
        void compute(const T& precision = std::numeric_limits<T>::epsilon());
        void compute_nonconvex(const T& precision = std::numeric_limits<T>::epsilon());
        /* Getters */
        [[nodiscard]] const VectorND<T>& getSolution() const noexcept { return x; }
        [[nodiscard]] T value() const;
    private:
        void updateVariables(const VectorND<T>& direction);
        [[nodiscard]] T nextStepSize(const VectorND<T>& direction, size_t& blockedAt);
        void updateActiveConstraints(ConstraintMatrix& activeConstraints);
    };

    template<Scalar T>
    QuadraticProgramming<T>::QuadraticProgramming(const Matrix auto& objectiveMatG_,
                                                  const Vector auto& objectiveVecC_,
                                                  const Matrix auto& equalityConstraint_,
                                                  const Matrix auto& inequalityConstraint_,
                                                  const Vector auto& initial)
            : objectiveMatG(objectiveMatG_)
            , objectiveVecC(objectiveVecC_)
            , x(initial)
            , activeConstraintFlags(Array<bool>(inequalityConstraint_.getRow(), false)) {
        assert(objectiveMatG.getRow() == objectiveVecC.getLength());
        assert(equalityConstraint_.getCol() == 0 || equalityConstraint_.getCol() == objectiveVecC.getLength() + 1);
        assert(inequalityConstraint_.getCol() == 0 || inequalityConstraint_.getCol() == objectiveVecC.getLength() + 1);
        if (equalityConstraint_.getSize() > 0)
            equalityConstraint = equalityConstraint_;
        if (inequalityConstraint_.getSize() > 0)
            inequalityConstraint = inequalityConstraint_;
    }

    template<Scalar T>
    void QuadraticProgramming<T>::compute(const T& precision) {
        ConstraintMatrix activeConstraints = equalityConstraint;
        while (true) {
            const EqualityQuadraticProgramming<T> EQP(objectiveMatG, objectiveVecC, activeConstraints, x);
            const VectorND<T> vec_p = EQP.getSolution() - x;
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
    template<Scalar T>
    void QuadraticProgramming<T>::compute_nonconvex(const T& precision) {
        ConstraintMatrix activeConstraints = equalityConstraint;
        while (true) {
            const EqualityQuadraticProgramming<T> EQP(objectiveMatG, objectiveVecC, activeConstraints, x);
            const VectorND<T> vec_p = EQP.getSolution() - x;
            if (vec_p.norm() <= x.norm() * precision)
                break;
            else
                updateVariables(vec_p);
            updateActiveConstraints(activeConstraints);
        };
    }

    template<Scalar T>
    T QuadraticProgramming<T>::value() const {
        return (x.copyToRowMatrix() * objectiveMatG).compute().row(0) * x * T(0.5) + objectiveVecC * x;
    }

    template<Scalar T>
    void QuadraticProgramming<T>::updateVariables(const VectorND<T>& direction) {
        size_t blockedAt = 0;
        const T step = nextStepSize(direction, blockedAt);
        x = x + step * direction;
        if (step != T(1)) {
            assert(activeConstraintFlags[blockedAt] == false);
            activeConstraintFlags[blockedAt] = true;
        }
    }

    template<Scalar T>
    T QuadraticProgramming<T>::nextStepSize(const VectorND<T>& direction, size_t& blockedAt) {
        T result = T(1);
        for (size_t i = 0; i < activeConstraintFlags.getLength(); ++i) {
            const bool isActive = activeConstraintFlags[i];
            if (!isActive) {
                const auto row = inequalityConstraint.row(i);
                const auto head = row.head(row.getLength() - 1);
                const T dot = head * direction;
                if (dot.isNegative()) {
                    const T stepSize_i = (row[row.getLength() - 1] - head * x) / dot;
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

    template<Scalar T>
    void QuadraticProgramming<T>::updateActiveConstraints(ConstraintMatrix& activeConstraints) {
        const size_t activeInequality = std::count(activeConstraintFlags.begin(), activeConstraintFlags.end(), true);
        activeConstraints.resize(equalityConstraint.getRow() + activeInequality, x.getLength() + 1);

        size_t activeInequalityIndex = 0;
        for (size_t i = 0; i < activeConstraintFlags.getLength(); ++i) {
            if (activeConstraintFlags[i]) {
                auto row = activeConstraints.row(equalityConstraint.getRow() + activeInequalityIndex);
                row = inequalityConstraint.row(i);
                ++activeInequalityIndex;
            }
        }
    }
}
