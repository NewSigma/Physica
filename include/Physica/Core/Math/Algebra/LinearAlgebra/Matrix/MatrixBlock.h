/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/LValueVector.h"

namespace Physica::Core {
    template<class MatrixType, size_t Row = Dynamic, size_t Column = Dynamic>
    class MatrixBlock;

    template<class MatrixType>
    class MatrixBlock<MatrixType, 1, Dynamic> : public LValueVector<MatrixBlock<MatrixType, 1, Dynamic>> {
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t row;
        size_t fromCol;
        size_t colCount;
    public:
        MatrixBlock(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : mat(mat_), row(row_), fromCol(fromCol_), colCount(colCount_) {}
        MatrixBlock(const MatrixBlock&) = delete;
        MatrixBlock(MatrixBlock&&) noexcept = delete;
        ~MatrixBlock() = default;
        /* Operators */
        MatrixBlock& operator=(const MatrixBlock&) = delete;
        MatrixBlock& operator=(MatrixBlock&&) noexcept = delete;
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < colCount); return mat(row, fromCol + index); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < colCount); return mat(row, fromCol + index); }
    };

    template<class MatrixType>
    class MatrixBlock<MatrixType, Dynamic, 1> : public LValueVector<MatrixBlock<MatrixType, Dynamic, 1>> {
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t col;
    public:
        MatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : mat(mat_), fromRow(fromRow_), rowCount(rowCount_), col(col_) {}
        MatrixBlock(const MatrixBlock&) = delete;
        MatrixBlock(MatrixBlock&&) noexcept = delete;
        ~MatrixBlock() = default;
        /* Operators */
        MatrixBlock& operator=(const MatrixBlock&) = delete;
        MatrixBlock& operator=(MatrixBlock&&) noexcept = delete;
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < rowCount); return mat(fromRow + index, col); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < rowCount); return mat(fromRow + index, col); }
    };

    template<class MatrixType>
    class MatrixBlock<MatrixType, Dynamic, Dynamic> {
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        MatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        MatrixBlock(const MatrixBlock&) = delete;
        MatrixBlock(MatrixBlock&&) noexcept = delete;
        ~MatrixBlock() = default;
        /* Operators */
        MatrixBlock& operator=(const MatrixBlock&) = delete;
        MatrixBlock& operator=(MatrixBlock&&) noexcept = delete;
        [[nodiscard]] ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t col) const;
    };

    template<class MatrixType>
    MatrixBlock<MatrixType, Dynamic, Dynamic>::MatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_)
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getColumn());
    }

    template<class MatrixType>
    typename MatrixBlock<MatrixType, Dynamic, Dynamic>::ScalarType&
    MatrixBlock<MatrixType, Dynamic, Dynamic>::operator()(size_t row, size_t col) {
        assert(row < rowCount);
        assert(col < colCount);
        return mat(row + fromRow, col + fromCol);
    }

    template<class MatrixType>
    const typename MatrixBlock<MatrixType, Dynamic, Dynamic>::ScalarType&
    MatrixBlock<MatrixType, Dynamic, Dynamic>::operator()(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat(row + fromRow, col + fromCol);
    }
}