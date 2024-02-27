/*
 * Copyright 2021-2023 WeiBo He.
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

namespace Physica::Core {
    //////////////////////////////////////////////Column-Element//////////////////////////////////////////////
    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>::rowSwap(size_t r1, size_t r2) {
        Derived& matrix = getDerived();
        const size_t row = matrix.getRow();
        const size_t column = matrix.getColumn();
        assert(r1 < row && r2 < row);
        for (size_t i = 0, temp = 0; i < column; ++i, temp += row)
            swap(matrix[temp + r1], matrix[temp + r2]);
    }

    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Element>::columnSwap(size_t c1, size_t c2) {
        Derived& matrix = getDerived();
        const size_t row = matrix.getRow();
        [[maybe_unused]] const size_t column = matrix.getColumn();
        assert(c1 < column && c2 < column);
        const size_t offset1 = c1 * row;
        const size_t offset2 = c2 * row;
        for (size_t i = 0; i < row; ++i)
            matrix[offset1 + i].swap(matrix[offset2 + i]);
    }
    //////////////////////////////////////////////Row-Element//////////////////////////////////////////////
    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>::rowSwap(size_t r1, size_t r2) {
        Derived& matrix = getDerived();
        [[maybe_unused]] const size_t row = matrix.getRow();
        const size_t column = matrix.getColumn();
        assert(r1 < row && r2 < row);
        const size_t offset1 = r1 * column;
        const size_t offset2 = r2 * column;
        for (size_t i = 0; i < column; ++i)
            matrix[offset1 + i].swap(matrix[offset2 + i]);
    }

    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Element>::columnSwap(size_t c1, size_t c2) {
        Derived& matrix = getDerived();
        const size_t row = matrix.getRow();
        const size_t column = matrix.getColumn();
        assert(c1 < column && c2 < column);
        for (size_t i = 0, temp = 0; i < column; ++i, temp += row)
            swap(matrix[temp + c1], matrix[temp + c2]);
    }
    //////////////////////////////////////////////Column-Vector//////////////////////////////////////////////
    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>::resize(size_t row, size_t column) {
        array.resize(column);
        for (auto& vector : array)
            vector.resize(row);
    }

    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>::rowSwap(size_t r1, size_t r2) {
        [[maybe_unused]] const size_t row = Base::getDerived().getRow();
        assert(r1 < row && r2 < row);
        for (auto& columnVector : array)
            columnVector[r1].swap(columnVector[r2]);
    }

    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Column | MatrixOption::Vector>::columnSwap(size_t c1, size_t c2) {
        [[maybe_unused]] const size_t column = Base::getDerived().getColumn();
        assert(c1 < column && c2 < column);
        array[c1].swap(array[c2]);
    }
    //////////////////////////////////////////////Row-Vector//////////////////////////////////////////////
    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>::resize(size_t row, size_t column) {
        array.resize(row);
        for (auto& vector : array)
            vector.resize(column);
    }

    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>::rowSwap(size_t r1, size_t r2) {
        [[maybe_unused]] const size_t row = Base::getDerived().getRow();
        assert(r1 < row && r2 < row);
        array[r1].swap(array[r2]);
    }

    template<class Derived>
    void DenseMatrixStorage<Derived, MatrixOption::Row | MatrixOption::Vector>::columnSwap(size_t c1, size_t c2) {
        [[maybe_unused]] const size_t column = Base::getDerived().getColumn();
        assert(c1 < column && c2 < column);
        for (auto& rowVector : array)
            rowVector[c1].swap(rowVector[c2]);
    }
}
