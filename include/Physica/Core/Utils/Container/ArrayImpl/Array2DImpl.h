/*
 * Copyright 2021-2025 Weibo He.
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

#include "../Array2D.h"

namespace Physica {
    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    Array2D<T, Option, Row, Col, Allocator>::Array2D(size_t row, size_t col) : r(row) {
        if constexpr (isVectorStorage) {
            if constexpr (isColMajor)
                arr = ArrayType(col, row);
            else
                arr = ArrayType(row, col);
        }
        else
            arr = ArrayType(row * col);
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    Array2D<T, Option, Row, Col, Allocator>::Array2D(size_t row, size_t col, T value) : r(row) {
        if constexpr (isVectorStorage) {
            if constexpr (isColMajor)
                arr = ArrayType(col, row, std::move(value));
            else
                arr = ArrayType(row, col, std::move(value));
        }
        else
            arr = ArrayType(row * col, std::move(value));
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    Array2D<T, Option, Row, Col, Allocator>::Array2D(std::initializer_list<InitializerType> list) : arr(list) {
        if constexpr (isVectorStorage) {
            const size_t length = arr.getLength();
            r = isColMajor ? (getSize() / length) : length;
        }
        else
            assert((Row != Dynamic || Col != Dynamic) && "[Error]: Either row or col must be given at compile");
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    Array2D<T, Option, Row, Col, Allocator>::Array2D(ArrayType arr_, IndexType r_) : arr(std::move(arr_)), r(r_) {}

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    T& Array2D<T, Option, Row, Col, Allocator>::operator()(size_t r, size_t c) {
        assert(r < getRow() && c < getCol());
        if constexpr (isVectorStorage) {
            if constexpr (isColMajor)
                return arr[c][r];
            else
                return arr[r][c];
        }
        else
            return arr[toIndex(r, c)];
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    const T& Array2D<T, Option, Row, Col, Allocator>::operator()(size_t r, size_t c) const {
        return const_cast<This&>(*this).operator()(r, c);
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    template<class... Args>
    void Array2D<T, Option, Row, Col, Allocator>::resize(size_t row, size_t col, Args&&... args) {
        assert((Row == row || Row == Dynamic) && "[Error]: Cannot resize a fixed array");
        assert((Col == col || Col == Dynamic) && "[Error]: Cannot resize a fixed array");
        if constexpr (isVectorStorage) {
            arr.resize(isColMajor ? col : row);
            for (auto& elem : arr)
                elem.resize(isColMajor ? row : col, std::forward<Args>(args)...);
        }
        else
            arr.resize(row * col, std::forward<Args>(args)...);
        r = row;
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    template<class... Args>
    void Array2D<T, Option, Row, Col, Allocator>::resize(size_t order) {
        resize(order, order);
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    void Array2D<T, Option, Row, Col, Allocator>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        arr.swap(obj.arr);
        std::swap(r, obj.r);
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    void Array2D<T, Option, Row, Col, Allocator>::swap_row(size_t r1, size_t r2) noexcept {
        assert(r1 < getRow() && r2 < getRow());
        if constexpr (isVectorStorage) {
            if constexpr (isColMajor)
                for (auto& elem : arr)
                    elem[r1].swap(elem[r2]);
            else
                arr[r1].swap(arr[r2]);
        }
        else {
            if constexpr (isColMajor) {
                const size_t row = getRow();
                const size_t col = getCol();
                for (size_t i = 0, temp = 0; i < col; ++i, temp += row)
                    arr[temp + r1].swap(arr[temp + r2]);
            }
            else {
                const size_t col = getCol();
                const size_t offset1 = r1 * col;
                const size_t offset2 = r2 * col;
                for (size_t i = 0; i < col; ++i)
                    arr[offset1 + i].swap(arr[offset2 + i]);
            }
        }
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    void Array2D<T, Option, Row, Col, Allocator>::swap_col(size_t c1, size_t c2) noexcept {
        assert(c1 < getCol() && c2 < getCol());
        if constexpr (isVectorStorage) {
            if constexpr (isColMajor) {
                arr[c1].swap(arr[c2]);
            }
            else
                for (auto& elem : arr)
                    elem[c1].swap(elem[c2]);
        }
        else {
            if constexpr (isColMajor) {
                const size_t row = getRow();
                const size_t offset1 = c1 * row;
                const size_t offset2 = c2 * row;
                for (size_t i = 0; i < row; ++i)
                    arr[offset1 + i].swap(arr[offset2 + i]);
            }
            else {
                const size_t row = getRow();
                const size_t col = getCol();
                for (size_t i = 0, temp = 0; i < col; ++i, temp += row)
                    arr[temp + c1].swap(arr[temp + c2]);
            }
        }
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    __host__ __device__ size_t Array2D<T, Option, Row, Col, Allocator>::getRow() const noexcept {
        if constexpr (Row == Dynamic) {
            if constexpr (Col == Dynamic)
                return r;
            else {
                const size_t size = getSize();
                assert(size % Col == 0);
                return size / Col;
            }
        }
        else
            return Row;
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    __host__ __device__ size_t Array2D<T, Option, Row, Col, Allocator>::getCol() const noexcept {
        if constexpr (Col == Dynamic) {
            const size_t size = getSize();
            if constexpr (Row == Dynamic) {
                assert(r == 0 || size % getRow() == 0);
                return r == 0 ? 0 : size / getRow();
            }
            else {
                assert(size % Row == 0);
                return size / Row;
            }
        }
        else
            return Col;
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    __host__ __device__ size_t Array2D<T, Option, Row, Col, Allocator>::getSize() const noexcept {
        if constexpr (isVectorStorage)
            return arr.getLength() == 0 ? 0 : arr.getLength() * arr[0].getLength();
        else
            return arr.getLength();
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    __host__ __device__ T* Array2D<T, Option, Row, Col, Allocator>::data_ptr(size_t row, size_t col) {
        if constexpr (isVectorStorage) {
            const size_t major = isColMajor ? col : row;
            const size_t minor = isColMajor ? row : col;
            return arr[major].data() + minor;
        }
        else
            return arr.data() + toIndex(row, col);
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    __host__ __device__ const T* Array2D<T, Option, Row, Col, Allocator>::data_ptr(size_t row, size_t col) const {
        return const_cast<This&>(*this).data_ptr(row, col);
    }

    template<class T, int Option, size_t Row, size_t Col, class Allocator>
    __host__ __device__ size_t Array2D<T, Option, Row, Col, Allocator>::toIndex(size_t r, size_t c) const {
        if constexpr (isColMajor)
            return getRow() * c + r;
        else
            return getCol() * r + c;
    }
}
