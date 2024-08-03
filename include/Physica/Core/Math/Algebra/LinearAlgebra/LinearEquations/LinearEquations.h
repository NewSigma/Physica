/*
 * Copyright 2020-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h>

namespace Physica::Core {
    /**
     * Unfinished:
     * If the equations do not have the unique solution, the program will throw a divide zero exception and stop.
     * Change the bias in LU method to solve a family of equations.
     *
     * References:
     * [1] William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery. C++数值算法(第二版)[M]. 北京: 电子工业出版社, 2005
     */
    template<class T, int Type = MatrixOption::Column | MatrixOption::Vector, size_t MaxRow = Dynamic, size_t MaxColumn = Dynamic>
    class LinearEquations {
        using This = LinearEquations<T, Type, MaxRow, MaxColumn>;

        DenseMatrix<T, Type, MaxRow, MaxColumn> working;
    public:
        explicit LinearEquations(DenseMatrix<T, Type, MaxRow, MaxColumn> working_);
        LinearEquations(const This& l) = default;
        LinearEquations(This&& l) noexcept = default;
        ~LinearEquations() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        void gaussJordanPartial();
        void gaussJordanComplete();
        void gaussEliminationPartial();
        void gaussEliminationComplete();
        void swap(This& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] const DenseMatrix<T, Type, MaxRow, MaxColumn>& getWorking() const noexcept { return working; }
        [[nodiscard]] auto getSolution() const noexcept { return working.col(working.getColumn() - 1); }
    private:
        void upperEliminate(size_t index);
        void lowerEliminate(size_t index);
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
