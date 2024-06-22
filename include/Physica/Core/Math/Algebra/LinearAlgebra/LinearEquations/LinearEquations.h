/*
 * Copyright 2020-2023 WeiBo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixOperation.h"

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
    template<class T, int Type = MatrixOption::Column | MatrixOption::Vector
            , size_t MaxRow = Utils::Dynamic, size_t MaxColumn = Utils::Dynamic>
    class LinearEquations {
        using Operation = MatrixOperation<T, Type, MaxRow, MaxColumn>;

        DenseMatrix<T, Type, MaxRow, MaxColumn> working;
    public:
        explicit LinearEquations(DenseMatrix<T, Type, MaxRow, MaxColumn> working_);
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
        void swap(LinearEquations& __restrict equ) noexcept;
        /* Getters */
        [[nodiscard]] const DenseMatrix<T, Type, MaxRow, MaxColumn>& getWorking() const noexcept { return working; }
        [[nodiscard]] auto getSolution() { return working.col(working.getColumn() - 1); }
    };
}

namespace std {
    template<class T, int Type, size_t MaxRow, size_t MaxColumn>
    inline void swap(Physica::Core::LinearEquations<T, Type, MaxRow, MaxColumn>& __restrict equ1,
                     Physica::Core::LinearEquations<T, Type, MaxRow, MaxColumn>& __restrict equ2) noexcept {
        equ1.swap(equ2);
    }
}

#include "LinearEquationsImpl.h"
