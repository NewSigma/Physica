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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Matrix/MatrixOperation.h"

namespace Physica::Core {
    /**
     * Unfinished:
     * If the equations do not have the unique solution, the program will throw a divide zero exception and stop.
     * Change the bias in LU method to solve a family of equations.
     * 
     * References:
     * [1] H.Press, William, A.Teukolsky, Saul, Vetterling, William T., Flannery, Brian P..
     * C++数值算法[M].北京: Publishing House of Electronics Industry, 2009
     */
    template<class T = MultiScalar, int type = DenseMatrixOption::Column | DenseMatrixOption::Vector
            , size_t maxRow = Utils::Dynamic, size_t maxColumn = Utils::Dynamic>
    class LinearEquations {
        using Operation = MatrixOperation<T, type, maxRow, maxColumn>;

        DenseMatrix<T, type, maxRow, maxColumn> working;
    public:
        explicit LinearEquations(DenseMatrix<T, type, maxRow, maxColumn> working_);
        LinearEquations(const LinearEquations& l) = default;
        LinearEquations(LinearEquations&& l) noexcept = default;
        ~LinearEquations() = default;
        /* Operators */
        LinearEquations& operator=(LinearEquations l) noexcept;
        /* Operations */
        void gaussJordanPartial();
        void gaussJordanComplete();
        void gaussEliminationPartial();
        void gaussEliminationComplete();
        void lu();
        void conjugateGradient();
        /* Getters */
        [[nodiscard]] const DenseMatrix<T, type, maxRow, maxColumn>& getWorking() const noexcept { return working; }
        [[nodiscard]] auto getSolution() { return working.col(working.getColumn() - 1); }
        /* Helpers */
        void swap(LinearEquations& equ) noexcept;
    };

    template<class T, int type, size_t maxRow, size_t maxColumn>
    inline void swap(LinearEquations<T, type, maxRow, maxColumn>& equ1, LinearEquations<T, type, maxRow, maxColumn>& equ2) noexcept {
        equ1.swap(equ2);
    }
}

#include "LinearEquationsImpl.h"
