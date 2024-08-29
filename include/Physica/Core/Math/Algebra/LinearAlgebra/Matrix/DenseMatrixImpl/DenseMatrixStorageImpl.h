/*
 * Copyright 2021-2024 Weibo He.
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
    #define tparams class T, size_t Row, size_t Column, size_t MaxRow, size_t MaxColumn, class Allocator
    #define ColumnElementStorage DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>

    template<tparams>
    template<class... Args>
    void ColumnElementStorage::resize(size_t row, size_t column, Args&&... args) {
        Base::resize(row * column, std::forward<Args>(args)...);
        Dim::resize(row, column);
    }

    template<tparams>
    void ColumnElementStorage::rowSwap(size_t r1, size_t r2) {
        const size_t row = getRow();
        const size_t column = getColumn();
        assert(r1 < row && r2 < row);
        for (size_t i = 0, temp = 0; i < column; ++i, temp += row)
            swap((*this)[temp + r1], (*this)[temp + r2]);
    }

    template<tparams>
    void ColumnElementStorage::columnSwap(size_t c1, size_t c2) {
        const size_t row = getRow();
        [[maybe_unused]] const size_t column = getColumn();
        assert(c1 < column && c2 < column);
        const size_t offset1 = c1 * row;
        const size_t offset2 = c2 * row;
        for (size_t i = 0; i < row; ++i)
            (*this)[offset1 + i].swap((*this)[offset2 + i]);
    }

    #undef ColumnElementStorage
    #define RowElementStorage DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Column, MaxRow, MaxColumn, Allocator>

    template<tparams>
    template<class... Args>
    void RowElementStorage::resize(size_t row, size_t column, Args&&... args) {
        Base::resize(row * column, std::forward<Args>(args)...);
        Dim::resize(row, column);
    }

    template<tparams>
    void RowElementStorage::rowSwap(size_t r1, size_t r2) {
        [[maybe_unused]] const size_t row = getRow();
        const size_t column = getColumn();
        assert(r1 < row && r2 < row);
        const size_t offset1 = r1 * column;
        const size_t offset2 = r2 * column;
        for (size_t i = 0; i < column; ++i)
            (*this)[offset1 + i].swap((*this)[offset2 + i]);
    }

    template<tparams>
    void RowElementStorage::columnSwap(size_t c1, size_t c2) {
        const size_t row = getRow();
        const size_t column = getColumn();
        assert(c1 < column && c2 < column);
        for (size_t i = 0, temp = 0; i < column; ++i, temp += row)
            swap((*this)[temp + c1], (*this)[temp + c2]);
    }

    #undef RowElementStorage
    #define ColumnVectorStorage DenseMatrixStorage<T, MatrixOption::Column | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>

    template<tparams>
    template<class... Args>
    void ColumnVectorStorage::resize(
            size_t row, size_t column, Args&&... args) {
        array.resize(column);
        for (auto& vector : array)
            vector.resize(row, std::forward<Args>(args)...);
        Dim::resize(row, column);
    }

    template<tparams>
    void ColumnVectorStorage::rowSwap(size_t r1, size_t r2) {
        assert(r1 < getRow() && r2 < getRow());
        for (auto& columnVector : array)
            columnVector[r1].swap(columnVector[r2]);
    }

    template<tparams>
    void ColumnVectorStorage::columnSwap(size_t c1, size_t c2) {
        assert(c1 < getColumn() && c2 < getColumn());
        array[c1].swap(array[c2]);
    }

    #undef ColumnVectorStorage
    #define RowVectorStorage DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Column, MaxRow, MaxColumn, Allocator>

    template<tparams>
    template<class... Args>
    void RowVectorStorage::resize(size_t row, size_t column, Args&&... args) {
        array.resize(row);
        for (auto& vector : array)
            vector.resize(column, std::forward<Args>(args)...);
        Dim::resize(row, column);
    }

    template<tparams>
    void RowVectorStorage::rowSwap(size_t r1, size_t r2) {
        assert(r1 < getRow() && r2 < getRow());
        array[r1].swap(array[r2]);
    }

    template<tparams>
    void RowVectorStorage::columnSwap(size_t c1, size_t c2) {
        assert(c1 < getColumn() && c2 < getColumn());
        for (auto& rowVector : array)
            rowVector[c1].swap(rowVector[c2]);
    }

    #undef RowVectorStorage
    #undef tparam
}
