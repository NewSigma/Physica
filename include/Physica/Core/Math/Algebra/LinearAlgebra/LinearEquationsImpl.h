/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_LINEAREQUATIONSIMPL_H
#define PHYSICA_LINEAREQUATIONSIMPL_H

#include "LUDecomposition.h"
#include "Matrix/MatrixOperation.h"

namespace Physica::Core {
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>::LinearEquations(const DenseMatrix<T, type, maxRow, maxColumn>& m)
            : matrix(m) {}

    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>::LinearEquations(DenseMatrix<T, type, maxRow, maxColumn>&& m) noexcept
            : matrix(std::move(m)) {}

    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>::LinearEquations(LinearEquations&& l) noexcept
            : matrix(std::move(l.matrix)) {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
    }

    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>& LinearEquations<T, type, maxRow, maxColumn>::operator=(
            LinearEquations&& l) noexcept {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)
        matrix = std::move(l.matrix);
    }
    /*!
    * After LinearEquations is constructed, a method should be selected and called.
    * Finally, call getResult() to acquire the solution.
    *
    * Requirement:
    * Matrix should have the rank of n * (n + 1).
    * Matrix must not be singular or a divide zero exception will be thrown.
    *
    * Unfinished:
    * If the equations do not have the unique solution, the program will throw a divide zero exception and stop.
    * Change the bias in LU method to solve a family of equations.
    */
    //!Reference: Numerical Recipes in C++
    template<class T, DenseMatrixType type, size_t maxRow, size_t maxColumn>
    void LinearEquations<T, type, maxRow, maxColumn>::solve(LinearEquationsMethod method) {
        Q_UNUSED(type)
        Q_UNUSED(maxRow)
        Q_UNUSED(maxColumn)

        typedef MatrixOperation<T, type, maxRow, maxColumn> Operation;
        const auto rank = matrix.getRow();
        Q_ASSERT(rank + 1 == matrix.getColumn());
        switch(method) {
            case GaussJordanPartial:
                for(size_t i = 0; i < rank; ++i) {
                    Operation::partialPivoting(matrix, i);
                    Operation::upperEliminate(matrix, i);
                    Operation::lowerEliminate(matrix, i);
                }
                for(size_t i = 0; i < rank; ++i)
                    matrix(i, rank) /= matrix(i, i);
                break;
            case GaussJordanComplete:
                for(size_t i = 0; i < rank; ++i) {
                    Operation::completePivoting(matrix, i);
                    Operation::upperEliminate(matrix, i);
                    Operation::lowerEliminate(matrix, i);
                }
                for(size_t i = 0; i < rank; ++i)
                    matrix(i, rank) /= matrix(i, i);
                break;
            case GaussEliminationPartial:
                for(size_t i = 0; i < rank; ++i) {
                    Operation::partialPivoting(matrix, i);
                    Operation::lowerEliminate(matrix, i);
                }
                for(size_t i = rank - 1; i > 0; --i) {
                    matrix(i, rank) /= matrix(i, i);
                    for(size_t j = 0; j < i; ++j)
                        matrix(j, rank) -= matrix(j, i) * matrix(i, rank);
                }
                matrix(0, rank) /= matrix(0, 0);
                break;
            case GaussEliminationComplete:
                for(size_t i = 0; i < rank; ++i) {
                    Operation::completePivoting(matrix, i);
                    Operation::lowerEliminate(matrix, i);
                }
                for(size_t i = rank - 1; i > 0; --i) {
                    matrix(i, rank) /= matrix(i, i);
                    for(size_t j = 0; j < i; ++j)
                        matrix(j, rank) -= matrix(j, i) * matrix(i, rank);
                }
                matrix(0, rank) /= matrix(0, 0);
                break;
            case LUMethod:
                LUDecomposition lu(std::move(matrix));
                matrix = lu.release();
                for(size_t i = 0; i < rank - 1; ++i) {
                    for(size_t j = i + 1; j < rank; ++j)
                        matrix(j, rank) -= matrix(j, i) * matrix(i, rank);
                }
                for(size_t i = rank - 1; i > 0; --i) {
                    matrix(i, rank) /= matrix(i, i);
                    for(size_t j = 0; j < i; ++j)
                        matrix(j, rank) -= matrix(j, i) * matrix(i, rank);
                }
                matrix(0, rank) /= matrix(0, 0);
                break;
        }
    }
}

#endif