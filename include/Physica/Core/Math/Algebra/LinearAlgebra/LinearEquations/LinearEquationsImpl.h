/*
 * Copyright 2020-2024 WeiBo He.
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
    template<class T, int Type, size_t MaxRow, size_t MaxColumn>
    LinearEquations<T, Type, MaxRow, MaxColumn>::LinearEquations(DenseMatrix<T, Type, MaxRow, MaxColumn> working_)
            : working(std::move(working_)) {
        assert(working.getRow() + 1 == working.getColumn());
    }

    template<class T, int Type, size_t MaxRow, size_t MaxColumn>
    void LinearEquations<T, Type, MaxRow, MaxColumn>::gaussJordanPartial() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            Operation::partialPivoting(working, i);
            Operation::upperEliminate(working, i);
            Operation::lowerEliminate(working, i);
        }
        for (size_t i = 0; i < rank; ++i)
            working(i, rank) /= working(i, i);
    }

    template<class T, int Type, size_t MaxRow, size_t MaxColumn>
    void LinearEquations<T, Type, MaxRow, MaxColumn>::gaussJordanComplete() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            Operation::completePivoting(working, i);
            Operation::upperEliminate(working, i);
            Operation::lowerEliminate(working, i);
        }
        for (size_t i = 0; i < rank; ++i)
            working(i, rank) /= working(i, i);
    }

    template<class T, int Type, size_t MaxRow, size_t MaxColumn>
    void LinearEquations<T, Type, MaxRow, MaxColumn>::gaussEliminationPartial() {
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

    template<class T, int Type, size_t MaxRow, size_t MaxColumn>
    void LinearEquations<T, Type, MaxRow, MaxColumn>::gaussEliminationComplete() {
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

    template<class T, int Type, size_t MaxRow, size_t MaxColumn>
    void LinearEquations<T, Type, MaxRow, MaxColumn>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
    }
}
