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
    template<Matrix M>
    class device_obj<RMatrixBlock<M, 1, Dynamic>> : public device_obj<RValueVector<RMatrixBlock<M, 1, Dynamic>>> {
        using host_obj = RMatrixBlock<M, 1, Dynamic>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<M>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(const device_obj<M>& mat_, size_t fromRow_, size_t fromCol_, size_t colCount_)
                : mat(asStruct(mat_)), fromRow(fromRow_), fromCol(fromCol_), colCount(colCount_) {
            assert(fromRow < mat_.getRow());
            assert(fromCol + colCount <= mat_.getCol());
        }
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t index) const {
            assert(index < colCount);
            return mat.getDerived().calc(fromRow, fromCol + index);
        }
        [[nodiscard]] __device__ Tv calc_value(size_t index) const {
            assert(index < colCount);
            return mat.getDerived().calc_value(fromRow, fromCol + index);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return colCount; }
    };

    template<Matrix M>
    class device_obj<RMatrixBlock<M, Dynamic, 1>> : public device_obj<RValueVector<RMatrixBlock<M, Dynamic, 1>>> {
        using host_obj = RMatrixBlock<M, Dynamic, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
    public:
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<M>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(const device_obj<M>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(asStruct(mat_)), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat_.getRow());
            assert(fromCol < mat_.getCol());
        }
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t index) const {
            assert(index < rowCount);
            return mat.getDerived().calc(fromRow + index, fromCol);
        }
        [[nodiscard]] __device__ Tv calc_value(size_t index) const {
            assert(index < rowCount);
            return mat.getDerived().calc_value(fromRow + index, fromCol);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return rowCount; }
    };

    template<Matrix M>
    class device_obj<RMatrixBlock<M, Dynamic, Dynamic>> : public device_obj<RValueMatrix<RMatrixBlock<M, Dynamic, Dynamic>>> {
        using host_obj = RMatrixBlock<M, Dynamic, Dynamic>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<RMatrixBlock<M, Dynamic, Dynamic>>>;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<const device_obj<M>> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(const device_obj<M>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return colCount; }
    };

    template<Matrix M>
    __host__ __device__ device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::device_obj(
            const device_obj<M>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(asStruct(mat_))
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat_.getRow());
        assert((fromCol + colCount) <= mat_.getCol());
    }

    template<Matrix M>
    __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::calc(size_t row, size_t col) const -> T {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.getDerived().calc(row + fromRow, col + fromCol);
    }

    template<Matrix M>
    __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::calc_value(size_t row, size_t col) const -> Tv {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.getDerived().calc_value(row + fromRow, col + fromCol);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<device_obj<RMatrixBlock<M, Row, Col>>> : public Traits<RMatrixBlock<M, Row, Col>> {};
}
