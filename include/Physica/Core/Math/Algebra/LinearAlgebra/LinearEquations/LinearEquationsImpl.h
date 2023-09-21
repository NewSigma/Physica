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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/PLUDecomposition.h"

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
    void LinearEquations<T, type, maxRow, maxColumn>::swap(LinearEquations& equ) noexcept {
        assert(this != &equ && "[Error]: Self swap is likely a bug");
        working.swap(equ.working);
    }
}
