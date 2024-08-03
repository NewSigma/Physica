/*
 * Copyright 2023-2024 Weibo He.
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

namespace Physica::Core {
    template<class MatrixType>
    class device_obj<RowRVector<MatrixType>> : public device_obj<RValueVector<RowRVector<MatrixType>>> {
    public:
        using Base = device_obj<RValueVector<RowRVector<MatrixType>>>;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        device_obj<MatrixType>& mat;
        size_t row;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<RValueMatrix<MatrixType>>& mat_, size_t row_, size_t fromCol_, size_t colCount_)
                : mat(mat_.getDerived()), row(row_), fromCol(fromCol_), colCount(colCount_) {
            assert(row < mat.getRow());
            assert(fromCol + colCount <= mat.getColumn());
        }
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { assert(index < colCount); return mat.calc(row, fromCol + index); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return colCount; }
    };

    template<class MatrixType>
    class device_obj<ColRVector<MatrixType>> : public device_obj<RValueVector<ColRVector<MatrixType>>> {
    public:
        using Base = device_obj<RValueVector<ColRVector<MatrixType>>>;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        device_obj<MatrixType>& mat;
        size_t col;
        size_t fromRow;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<RValueMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : mat(mat_.getDerived()), col(col_), fromRow(fromRow_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(col < mat.getColumn());
        }
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t index) const { assert(index < rowCount); return mat.calc(fromRow + index, col); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return rowCount; }
    };

    template<class MatrixType>
    class device_obj<RMatrixBlock<MatrixType, 1, Dynamic>>
            : public device_obj<RValueMatrix<RMatrixBlock<MatrixType, 1, Dynamic>>>
            , public device_obj<RowRVector<MatrixType>> {
    public:
        using Base = device_obj<RValueMatrix<RMatrixBlock<MatrixType, 1, Dynamic>>>;
        using VectorBase = device_obj<RowRVector<MatrixType>>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        __host__ __device__ device_obj(device_obj<RValueMatrix<MatrixType>>& mat_, size_t row_, size_t fromCol_, size_t colCount_)
                : VectorBase(mat_, row_, fromCol_, colCount_) {}
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return VectorBase::calc(col); }
        using VectorBase::calc;
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return VectorBase::getLength(); }
        //using VectorBase::conjugate;  // Not implemented
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] __host__ __device__ Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class device_obj<RMatrixBlock<MatrixType, Dynamic, 1>>
            : public device_obj<RValueMatrix<RMatrixBlock<MatrixType, Dynamic, 1>>>
            , public device_obj<ColRVector<MatrixType>> {
    public:
        using Base = device_obj<RValueMatrix<RMatrixBlock<MatrixType, Dynamic, 1>>>;
        using VectorBase = device_obj<ColRVector<MatrixType>>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        __host__ __device__ device_obj(device_obj<RValueMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, [[maybe_unused]] size_t col) const { assert(col == 0); return VectorBase::calc(row); }
        using VectorBase::calc;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return VectorBase::getLength(); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getColumn() noexcept { return 1; }
        //using VectorBase::conjugate;  // Not implemented
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] __host__ __device__ Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class device_obj<RMatrixBlock<MatrixType, Dynamic, Dynamic>> : public device_obj<RValueMatrix<RMatrixBlock<MatrixType, Dynamic, Dynamic>>> {
    public:
        using Base = device_obj<RValueMatrix<RMatrixBlock<MatrixType, Dynamic, Dynamic>>>;
        using typename Base::ScalarType;
    private:
        device_obj<MatrixType>& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<RValueMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return colCount; }
    };

    template<class MatrixType>
    __host__ __device__ device_obj<RMatrixBlock<MatrixType, Dynamic, Dynamic>>::device_obj(
            device_obj<RValueMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_.getDerived())
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getColumn());
    }

    template<class MatrixType>
    __device__ typename device_obj<RMatrixBlock<MatrixType, Dynamic, Dynamic>>::ScalarType
    device_obj<RMatrixBlock<MatrixType, Dynamic, Dynamic>>::calc(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.calc(row + fromRow, col + fromCol);
    }
}

namespace Physica {
    using namespace Core;

    template<class MatrixType>
    class Traits<Core::device_obj<RowRVector<MatrixType>>> : public Traits<RowRVector<MatrixType>> {};

    template<class MatrixType>
    class Traits<Core::device_obj<ColRVector<MatrixType>>> : public Traits<ColRVector<MatrixType>> {};

    template<class MatrixType, size_t Row, size_t Column>
    class Traits<Core::device_obj<RMatrixBlock<MatrixType, Row, Column>>> : public Traits<RMatrixBlock<MatrixType, Row, Column>> {};
}
