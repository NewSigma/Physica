/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/SquareMatrix.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/LUDecomposition.h>
#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearEquations.h"

namespace Physica::Core {
    /*!
     * After LinearEquations is constructed, a method should be selected and called.
     * Finally, call getResult() to acquire the solution.
     *
     * Requirement:
     * Matrix should have the rank of n * (n + 1).
     * Matrix must not be singular or a divide zero exception will be thrown.
     *
     * Optimize:
     * If the equations do not have the unique solution, the program will throw a divide zero exception and stop.
     *
     * Unfinished:
     * Change the bias in LU method to solve a family of equations.
     */
    LinearEquations::LinearEquations(Matrix& m) noexcept : matrix(m) {}
    //!Reference: Numerical Recipes in C++
    const Vector<MultiScalar>& LinearEquations::solve(LinearEquationsMethod method) {
        const auto rank = matrix.row();
        Q_ASSERT(rank + 1 == matrix.column());
        switch(method) {
            case GaussJordanPartial:
                for(size_t i = 0; i < rank; ++i) {
                    matrix.partialPivoting(i);
                    matrix.upperEliminate(i);
                    matrix.lowerEliminate(i);
                }
                for(size_t i = 0; i < rank; ++i)
                    matrix(i, rank) /= matrix(i, i);
                break;
            case GaussJordanComplete:
                for(size_t i = 0; i < rank; ++i) {
                    matrix.completePivoting(i);
                    matrix.upperEliminate(i);
                    matrix.lowerEliminate(i);
                }
                for(size_t i = 0; i < rank; ++i)
                    matrix(i, rank) /= matrix(i, i);
                break;
            case GaussEliminationPartial:
                for(size_t i = 0; i < rank; ++i) {
                    matrix.partialPivoting(i);
                    matrix.lowerEliminate(i);
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
                    matrix.completePivoting(i);
                    matrix.lowerEliminate(i);
                }
                for(size_t i = rank - 1; i > 0; --i) {
                    matrix(i, rank) /= matrix(i, i);
                    for(size_t j = 0; j < i; ++j)
                        matrix(j, rank) -= matrix(j, i) * matrix(i, rank);
                }
                matrix(0, rank) /= matrix(0, 0);
                break;
            case LUMethod:
                SquareMatrix* square = matrix.getType() == Matrix::Column
                        ? static_cast<SquareMatrix*>(new ColumnSquareMatrix(std::move(matrix)))
                        : static_cast<SquareMatrix*>(new RowSquareMatrix(std::move(matrix)));
                LUDecomposition lu(*square);
                static_cast<CStyleArray<Vector<MultiScalar>, Dynamic, Dynamic>&>(matrix) = std::move(*square);
                delete square;
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