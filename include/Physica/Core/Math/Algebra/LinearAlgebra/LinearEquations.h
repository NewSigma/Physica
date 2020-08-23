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
#ifndef PHYSICA_LINEAREQUATIONS_H
#define PHYSICA_LINEAREQUATIONS_H

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Matrix.h"

namespace Physica::Core {
    /*!
     * Declare common parts for template instances of LinearEquations.
     */
    class AbstractLinearEquations {
    public:
        enum LinearEquationsMethod {
            GaussJordanPartial,
            GaussJordanComplete,
            GaussEliminationPartial,
            GaussEliminationComplete,
            LUMethod
        };
    };
    /*!
     * Solve linear equations.
     */
    template<class T = MultiScalar, MatrixType type = Column, size_t maxRow = Dynamic, size_t maxColumn = Dynamic>
    class LinearEquations : AbstractLinearEquations{
        Matrix<T, type, maxRow, maxColumn> matrix;
    public:
        explicit LinearEquations(const Matrix<T, type, maxRow, maxColumn>& m);
        explicit LinearEquations(Matrix<T, type, maxRow, maxColumn>&& m) noexcept;
        LinearEquations(const LinearEquations& l) = default;
        LinearEquations(LinearEquations&& l) noexcept;
        ~LinearEquations() = default;
        /* Operators */
        LinearEquations& operator=(const LinearEquations& l) = default;
        LinearEquations& operator=(LinearEquations&& l) noexcept;
        /* Operations */
        void solve(LinearEquationsMethod method);
        /* Helpers */
        [[nodiscard]] Matrix<T, type, maxRow, maxColumn>&& release() noexcept { return std::move(matrix); }
        /* Getters */
        [[nodiscard]] const Matrix<T, type, maxRow, maxColumn>& getMatrix() const noexcept { return matrix; }
        [[nodiscard]] T& getResult(size_t index) { return matrix(index, matrix.getColumn() - 1); }
    };
}

#include "LinearEquationsImpl.h"

#endif