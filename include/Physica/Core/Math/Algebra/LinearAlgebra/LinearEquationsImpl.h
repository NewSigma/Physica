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

#include <iostream>
#include "Matrix/MatrixDecomposition/PLUDecomposition.h"

namespace Physica::Core {
    template<class T, int type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>::LinearEquations(DenseMatrix<T, type, maxRow, maxColumn> working_)
            : working(std::move(working_)) {
        assert(working.getRow() + 1 == working.getColumn());
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>&
    LinearEquations<T, type, maxRow, maxColumn>::operator=(LinearEquations l) noexcept {
        swap(l);
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    void LinearEquations<T, type, maxRow, maxColumn>::gaussJordanPartial() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            Operation::partialPivoting(working, i);
            Operation::upperEliminate(working, i);
            Operation::lowerEliminate(working, i);
        }
        for (size_t i = 0; i < rank; ++i)
            working(i, rank) /= working(i, i);
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    void LinearEquations<T, type, maxRow, maxColumn>::gaussJordanComplete() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            Operation::completePivoting(working, i);
            Operation::upperEliminate(working, i);
            Operation::lowerEliminate(working, i);
        }
        for (size_t i = 0; i < rank; ++i)
            working(i, rank) /= working(i, i);
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    void LinearEquations<T, type, maxRow, maxColumn>::gaussEliminationPartial() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            Operation::partialPivoting(working, i);
            Operation::lowerEliminate(working, i);
        }
        for (size_t i = rank - 1; i > 0; --i) {
            working(i, rank) /= working(i, i);
            for (size_t j = 0; j < i; ++j)
                working(j, rank) -= working(j, i) * working(i, rank);
        }
        working(0, rank) /= working(0, 0);
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    void LinearEquations<T, type, maxRow, maxColumn>::gaussEliminationComplete() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            Operation::completePivoting(working, i);
            Operation::lowerEliminate(working, i);
        }
        for (size_t i = rank - 1; i > 0; --i) {
            working(i, rank) /= working(i, i);
            for (size_t j = 0; j < i; ++j)
                working(j, rank) -= working(j, i) * working(i, rank);
        }
        working(0, rank) /= working(0, 0);
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    void LinearEquations<T, type, maxRow, maxColumn>::lu() {
        const auto rank = working.getRow();
        PLUDecomposition lu(std::move(working));
        working = lu.release();
        for (size_t i = 0; i < rank - 1; ++i) {
            auto bottom = working.bottomRows(i + 1);
            auto col = bottom.col(rank);
            col -= bottom.col(i).asVector() * working(i, rank);
        }
        for (size_t i = rank - 1; i > 0; --i) {
            working(i, rank) /= working(i, i);
            auto top = working.topRows(i);
            auto col = top.col(rank);
            col -= top.col(i).asVector() * working(i, rank);
        }
        working(0, rank) /= working(0, 0);
    }
    /**
     * Pertinent for large-scale problem, coefficient matrix must be symmetric and positive definite
     * 
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006.112
     */
    template<class T, int type, size_t maxRow, size_t maxColumn>
    void LinearEquations<T, type, maxRow, maxColumn>::conjugateGradient() {
        using VectorType = Vector<T, maxRow>;
        
        const size_t row = working.getRow();
        const size_t column = working.getColumn();
        const auto mat_A = working.leftCols(column - 1);
        VectorType x = VectorType::Zeros(row);
        VectorType residual = -working.col(column - 1).asVector();
        VectorType p = working.col(column - 1);

        T squaredNorm = residual.squaredNorm();
        while(squaredNorm > std::numeric_limits<T>::epsilon()) {
            const VectorType temp = mat_A * p;
            const T step = squaredNorm / (p * temp);
            x += step * p;
            residual += step * temp;
            const T next_squaredNorm = residual.squaredNorm();
            const T beta = next_squaredNorm / squaredNorm;
            squaredNorm = next_squaredNorm;
            p = beta * p - residual;
        }
        auto result = working.col(column - 1);
        result = x;
    }

    template<class T, int type, size_t maxRow, size_t maxColumn>
    void LinearEquations<T, type, maxRow, maxColumn>::swap(LinearEquations& equ) noexcept {
        working.swap(equ.working);
    }
}
