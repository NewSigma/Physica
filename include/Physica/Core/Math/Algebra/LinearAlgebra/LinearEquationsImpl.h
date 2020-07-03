/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_LINEAREQUATIONSIMPL_H
#define PHYSICA_LINEAREQUATIONSIMPL_H

#include "LUDecomposition.h"
#include "Matrix/MatrixOperation.h"

namespace Physica::Core {
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>::LinearEquations(const Matrix<T, type, maxRow, maxColumn>& m)
            : matrix(m) {}

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>::LinearEquations(Matrix<T, type, maxRow, maxColumn>&& m) noexcept
            : matrix(std::move(m)) {}

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>::LinearEquations(LinearEquations&& l) noexcept
            : matrix(std::move(l.matrix)) {}

    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    LinearEquations<T, type, maxRow, maxColumn>& LinearEquations<T, type, maxRow, maxColumn>::operator=(
            LinearEquations&& l) noexcept {
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
    template<class T, MatrixType type, size_t maxRow, size_t maxColumn>
    const typename Matrix<T, type, maxRow, maxColumn>::VectorType&
    LinearEquations<T, type, maxRow, maxColumn>::solve(LinearEquationsMethod method) {
        typedef MatrixOperation<T, type, maxRow, maxColumn> Operation;
        const auto rank = matrix.row();
        Q_ASSERT(rank + 1 == matrix.column());
        switch(method) {
            case GaussJordanPartial:
                for(size_t i = 0; i < rank; ++i) {
                    Operation::partialPivoting(i);
                    Operation::upperEliminate(i);
                    Operation::lowerEliminate(i);
                }
                for(size_t i = 0; i < rank; ++i)
                    matrix(i, rank) /= matrix(i, i);
                break;
            case GaussJordanComplete:
                for(size_t i = 0; i < rank; ++i) {
                    Operation::completePivoting(i);
                    Operation::upperEliminate(i);
                    Operation::lowerEliminate(i);
                }
                for(size_t i = 0; i < rank; ++i)
                    matrix(i, rank) /= matrix(i, i);
                break;
            case GaussEliminationPartial:
                for(size_t i = 0; i < rank; ++i) {
                    Operation::partialPivoting(i);
                    Operation::lowerEliminate(i);
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
                    Operation::completePivoting(i);
                    Operation::lowerEliminate(i);
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
        return matrix[matrix.row()];
    }
}

#endif