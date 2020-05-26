/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include <Core/Header/SquareMatrix.h>
#include "Core/Header/LinearEquations.h"

namespace Physica::Core {
    /*
     * After LinearEquations is constructed, a method should be selected and called.
     * Finally, call getResult() to acquire the solution.
     *
     * Requirement:
     * Matrix should have the rank of n * (n + 1).
     * Matrix must not be singular or a divide zero exception will be thrown.
     *
     * Optimize:
     * If the equations do not have the unique solution, the program will throw a divide zero exception and stop.
     */
    LinearEquations::LinearEquations(Matrix*& m) noexcept : matrix(*m) {
        m = nullptr;
    }

    LinearEquations::LinearEquations(Matrix*&& m) noexcept : matrix(*m) {}

    LinearEquations::~LinearEquations() {
        delete &matrix;
    }

    const Vector& LinearEquations::solve(LinearEquationsMethod method) {
        const auto rank = matrix.row();
        Q_ASSERT(rank + 1 == matrix.column());
        switch(method) {
            case GaussJordanPartial:
                for(size_t i = 0; i < rank; ++i) {
                    matrix.partialPivoting(i);
                    matrix.upperEliminate(i);
                    matrix.lowerEliminate(i);
                }
                break;
            case GaussJordanComplete:
                for(size_t i = 0; i < rank; ++i) {
                    matrix.completePivoting(i);
                    matrix.upperEliminate(i);
                    matrix.lowerEliminate(i);
                }
                break;
            case GaussEliminationPartial:
                for(size_t i = 0; i < rank; ++i) {
                    matrix.partialPivoting(i);
                    matrix.lowerEliminate(i);
                }
                for(size_t i = rank - 1; i >= 0; --i)
                    for(size_t j = 0; j < i; ++j)
                        matrix(j, rank) -= matrix(j, i) * matrix(i, rank);
                break;
            case GaussEliminationComplete:
                for(size_t i = 0; i < rank; ++i) {
                    matrix.partialPivoting(i);
                    matrix.lowerEliminate(i);
                }
                for(size_t i = rank - 1; i >= 0; --i)
                    for(size_t j = 0; j < i; ++j)
                        matrix(j, rank) -= matrix(j, i) * matrix(i, rank);
                break;
            case LUDecomposition:
                for(size_t i = 0; i < rank; ++i) {
                    matrix.partialPivoting(i);
                    reinterpret_cast<SquareMatrix&>(matrix).LUDecompositionColumn(i); //NOLINT DataStructure not changed.
                }
                for(size_t i = rank - 1; i >= 0; --i)
                    for(size_t j = 0; j < i; ++j)
                        matrix(j, rank) -= matrix(j, i) * matrix(i, rank);
                for(size_t i = 0; i < rank; ++i)
                    for(size_t j = 0; j < i; ++j)
                        matrix(j, rank) -= matrix(j, i) * matrix(i, rank);
                break;
        }
        return matrix[matrix.row()];
    }
}