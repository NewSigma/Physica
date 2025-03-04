/*
 * Copyright 2023-2025 Weibo He.
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

#include "RMatrixBlock.h"

namespace Physica {
    template<Matrix T>
    class device_obj<RMatrixBlock<T, 1, Dynamic>> : public device_obj<RValueVector<RMatrixBlock<T, 1, Dynamic>>> {
        using host_obj = RMatrixBlock<T, 1, Dynamic>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    private:
        device_obj<T>& mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<T>& mat_, size_t fromRow_, size_t fromCol_, size_t colCount_)
                : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), colCount(colCount_) {
            assert(fromRow < mat.getRow());
            assert(fromCol + colCount <= mat.getCol());
        }
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            assert(index < colCount);
            return mat.calc(fromRow, fromCol + index);
        }
        [[nodiscard]] __device__ ValueType calc_value(size_t index) const {
            assert(index < colCount);
            return mat.calc_value(fromRow, fromCol + index);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return colCount; }
    };

    template<Matrix T>
    class device_obj<RMatrixBlock<T, Dynamic, 1>> : public device_obj<RValueVector<RMatrixBlock<T, Dynamic, 1>>> {
        using host_obj = RMatrixBlock<T, Dynamic, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    private:
        device_obj<T>& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<T>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const {
            assert(index < rowCount);
            return mat.calc(fromRow + index, fromCol);
        }
        [[nodiscard]] __device__ ValueType calc_value(size_t index) const {
            assert(index < rowCount);
            return mat.calc_value(fromRow + index, fromCol);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return rowCount; }
    };

    template<Matrix T>
    class device_obj<RMatrixBlock<T, Dynamic, Dynamic>> : public device_obj<RValueMatrix<RMatrixBlock<T, Dynamic, Dynamic>>> {
        using host_obj = RMatrixBlock<T, Dynamic, Dynamic>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<RMatrixBlock<T, Dynamic, Dynamic>>>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
    private:
        device_obj<T>& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<T>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ ValueType calc_value(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return colCount; }
    };

    template<Matrix T>
    __host__ __device__ device_obj<RMatrixBlock<T, Dynamic, Dynamic>>::device_obj(
            device_obj<T>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_)
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<Matrix T>
    __device__ auto device_obj<RMatrixBlock<T, Dynamic, Dynamic>>::calc(size_t row, size_t col) const -> ScalarType {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.calc(row + fromRow, col + fromCol);
    }

    template<Matrix T>
    __device__ auto device_obj<RMatrixBlock<T, Dynamic, Dynamic>>::calc_value(size_t row, size_t col) const -> ValueType {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.calc_value(row + fromRow, col + fromCol);
    }
}

namespace Physica {
    template<Matrix T, size_t Row, size_t Col>
    class Traits<device_obj<RMatrixBlock<T, Row, Col>>> : public Traits<RMatrixBlock<T, Row, Col>> {};
}
