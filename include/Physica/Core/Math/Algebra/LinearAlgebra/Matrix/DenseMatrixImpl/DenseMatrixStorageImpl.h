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
#define tparams class T, size_t Row, size_t Col, class Allocator
#define ColumnElementStorage DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Element, Row, Col, Allocator>

    template<tparams>
    template<class... Args>
    void ColumnElementStorage::resize(size_t row, size_t col, Args&&... args) {
        arr.resize(row * col, std::forward<Args>(args)...);
        Dim::resize(row, col);
    }

    template<tparams>
    void ColumnElementStorage::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Dim::swap(obj);
        arr.swap(obj.arr);
    }

    template<tparams>
    void ColumnElementStorage::swap_row(size_t r1, size_t r2) noexcept {
        const size_t row = getRow();
        const size_t col = getCol();
        assert(r1 < row && r2 < row);
        for (size_t i = 0, temp = 0; i < col; ++i, temp += row)
            arr[temp + r1].swap(arr[temp + r2]);
    }

    template<tparams>
    void ColumnElementStorage::swap_col(size_t c1, size_t c2) noexcept {
        const size_t row = getRow();
        [[maybe_unused]] const size_t col = getCol();
        assert(c1 < col && c2 < col);
        const size_t offset1 = c1 * row;
        const size_t offset2 = c2 * row;
        for (size_t i = 0; i < row; ++i)
            arr[offset1 + i].swap(arr[offset2 + i]);
    }

#undef ColumnElementStorage
#define RowElementStorage DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Element, Row, Col, Allocator>

    template<tparams>
    template<class... Args>
    void RowElementStorage::resize(size_t row, size_t col, Args&&... args) {
        arr.resize(row * col, std::forward<Args>(args)...);
        Dim::resize(row, col);
    }

    template<tparams>
    void RowElementStorage::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        Dim::swap(obj);
        arr.swap(obj.arr);
    }

    template<tparams>
    void RowElementStorage::swap_row(size_t r1, size_t r2) noexcept {
        [[maybe_unused]] const size_t row = getRow();
        const size_t col = getCol();
        assert(r1 < row && r2 < row);
        const size_t offset1 = r1 * col;
        const size_t offset2 = r2 * col;
        for (size_t i = 0; i < col; ++i)
            arr[offset1 + i].swap(arr[offset2 + i]);
    }

    template<tparams>
    void RowElementStorage::swap_col(size_t c1, size_t c2) noexcept {
        const size_t row = getRow();
        const size_t col = getCol();
        assert(c1 < col && c2 < col);
        for (size_t i = 0, temp = 0; i < col; ++i, temp += row)
            arr[temp + c1].swap(arr[temp + c2]);
    }

#undef RowElementStorage
#define ColumnVectorStorage DenseMatrixStorage<T, MatrixOption::Col | MatrixOption::Vector, Row, Col, Allocator>

    template<tparams>
    template<class... Args>
    void ColumnVectorStorage::resize(
            size_t row, size_t col, Args&&... args) {
        array.resize(col);
        for (auto& vector : array)
            vector.resize(row, std::forward<Args>(args)...);
        Dim::resize(row, col);
    }

    template<tparams>
    void ColumnVectorStorage::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        array.swap(obj.array);
        Dim::swap(obj);
    }

    template<tparams>
    void ColumnVectorStorage::swap_row(size_t r1, size_t r2) noexcept {
        assert(r1 < getRow() && r2 < getRow());
        for (auto& colVector : array)
            colVector[r1].swap(colVector[r2]);
    }

    template<tparams>
    void ColumnVectorStorage::swap_col(size_t c1, size_t c2) noexcept {
        assert(c1 < getCol() && c2 < getCol());
        array[c1].swap(array[c2]);
    }

#undef ColumnVectorStorage
#define RowVectorStorage DenseMatrixStorage<T, MatrixOption::Row | MatrixOption::Vector, Row, Col, Allocator>

    template<tparams>
    template<class... Args>
    void RowVectorStorage::resize(size_t row, size_t col, Args&&... args) {
        array.resize(row);
        for (auto& vector : array)
            vector.resize(col, std::forward<Args>(args)...);
        Dim::resize(row, col);
    }

    template<tparams>
    void RowVectorStorage::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        array.swap(obj.array);
        Dim::swap(obj);
    }

    template<tparams>
    void RowVectorStorage::swap_row(size_t r1, size_t r2) noexcept {
        assert(r1 < getRow() && r2 < getRow());
        array[r1].swap(array[r2]);
    }

    template<tparams>
    void RowVectorStorage::swap_col(size_t c1, size_t c2) noexcept {
        assert(c1 < getCol() && c2 < getCol());
        for (auto& rowVector : array)
            rowVector[c1].swap(rowVector[c2]);
    }

#undef RowVectorStorage
#undef tparams
}
