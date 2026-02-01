/*
 * Copyright 2020-2026 Weibo He.
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

namespace Physica {
    template<Scalar T>
    LinearSystem<T>::LinearSystem(WorkingMatrix&& working_) : working(std::move(working_)) {
        assert(working.getRow() + 1 == working.getCol() && "[Error]: Invalid working matrix");
    }

    template<Scalar T>
    LinearSystem<T>::LinearSystem(const Matrix auto& working_) : working(working_) {
        assert(working.getRow() + 1 == working.getCol() && "[Error]: Invalid working matrix");
    }

    template<Scalar T>
    void LinearSystem<T>::gaussJordanPartial() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            size_t pivot = working.leftCols(rank).pivotPartial(i);
            if (i != pivot)
                working.swap_row(i, pivot);

            upperEliminate(i);
            lowerEliminate(i);
        }
        for (size_t i = 0; i < rank; ++i)
            working[i, rank] /= working[i, i];
    }

    template<Scalar T>
    void LinearSystem<T>::gaussJordanComplete() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            size_t pivot = working.leftCols(rank).pivotComplete(i);
            if (i != pivot)
                working.swap_row(i, pivot);

            upperEliminate(i);
            lowerEliminate(i);
        }
        for (size_t i = 0; i < rank; ++i)
            working[i, rank] /= working[i, i];
    }

    template<Scalar T>
    void LinearSystem<T>::gaussEliminationPartial() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            size_t pivot = working.leftCols(rank).pivotPartial(i);
            if (i != pivot)
                working.swap_row(i, pivot);

            lowerEliminate(i);
        }
        for (size_t i = rank - 1; i > 0; --i) {
            working[i, rank] /= working[i, i];
            for (size_t j = 0; j < i; ++j)
                working[j, rank] -= working[j, i] * working[i, rank];
        }
        working[0, rank] /= working[0, 0];
    }

    template<Scalar T>
    void LinearSystem<T>::gaussEliminationComplete() {
        const auto rank = working.getRow();
        for (size_t i = 0; i < rank; ++i) {
            size_t pivot = working.leftCols(rank).pivotComplete(i);
            if (i != pivot)
                working.swap_row(i, pivot);

            lowerEliminate(i);
        }
        for (size_t i = rank - 1; i > 0; --i) {
            working[i, rank] /= working[i, i];
            for (size_t j = 0; j < i; ++j)
                working[j, rank] -= working[j, i] * working[i, rank];
        }
        working[0, rank] /= working[0, 0];
    }

    template<Scalar T>
    void LinearSystem<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        working.swap(obj.working);
    }

    template<Scalar T>
    void LinearSystem<T>::upperEliminate(size_t index) {
        assert(!working.operator[](index, index).isZero() && "[Error]: Matrix is singular");
        for(size_t i = 0; i < index; ++i)
            working.rowReduce(index, i, index);
    }

    template<Scalar T>
    void LinearSystem<T>::lowerEliminate(size_t index) {
        assert((!working[index, index].isZero()) && "[Error]: Matrix is singular");
        const auto r = working.getRow();
        for(size_t i = index + 1; i < r; ++i)
            working.rowReduce(index, i, index);
    }
}
