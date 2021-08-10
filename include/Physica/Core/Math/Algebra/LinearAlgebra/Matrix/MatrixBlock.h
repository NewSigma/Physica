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

    namespace Internal {
        template<class MatrixType, size_t Row, size_t Column>
        class Traits<MatrixBlock<MatrixType, Row, Column>> {
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static int MatrixOption = MatrixType::MatrixOption;
            constexpr static size_t RowAtCompile = Row;
            constexpr static size_t ColumnAtCompile = Column;
            constexpr static size_t MaxRowAtCompile = Row;
            constexpr static size_t MaxColumnAtCompile = Column;
            constexpr static size_t SizeAtCompile = Row * Column;
            constexpr static size_t MaxSizeAtCompile = SizeAtCompile;
        };
    }

    template<class MatrixType>
    class MatrixBlock<MatrixType, 1, Dynamic> : public LValueVector<MatrixBlock<MatrixType, 1, Dynamic>> {
    public:
        using Base = LValueVector<MatrixBlock<MatrixType, 1, Dynamic>>;
        using typename Base::ScalarType;
    private:
        MatrixType& mat;
        size_t row;
        size_t fromCol;
        size_t colCount;
    public:
        MatrixBlock(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : mat(mat_), row(row_), fromCol(fromCol_), colCount(colCount_) {}
        MatrixBlock(const MatrixBlock&) = default;
        MatrixBlock(MatrixBlock&&) noexcept = delete;
        ~MatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        MatrixBlock& operator=(const MatrixBlock& m) { Base::operator=(m); return *this; }
        MatrixBlock& operator=(MatrixBlock&& m) noexcept { Base::operator=(m); return *this; }
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < colCount); return mat(row, fromCol + index); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < colCount); return mat(row, fromCol + index); }
        /* Operations */
        void resize(size_t row, size_t col) { assert(row == 1 && col == colCount); }
        /* Getters */
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getColumn() const noexcept { return colCount; }
        [[nodiscard]] size_t getLength() const noexcept { return getColumn(); }
    };

    template<class MatrixType>
    class MatrixBlock<MatrixType, Dynamic, 1> : public LValueVector<MatrixBlock<MatrixType, Dynamic, 1>> {
    public:
        using Base = LValueVector<MatrixBlock<MatrixType, Dynamic, 1>>;
        using typename Base::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t col;
    public:
        MatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : mat(mat_), fromRow(fromRow_), rowCount(rowCount_), col(col_) {}
        MatrixBlock(const MatrixBlock&) = default;
        MatrixBlock(MatrixBlock&&) noexcept = delete;
        ~MatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        MatrixBlock& operator=(const MatrixBlock& m) { Base::operator=(m); return *this; }
        MatrixBlock& operator=(MatrixBlock&& m) noexcept { Base::operator=(m); return *this; }
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < rowCount); return mat(fromRow + index, col); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < rowCount); return mat(fromRow + index, col); }
        /* Operations */
        void resize(size_t row, size_t col) { assert(row == rowCount && col == 1); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return 1; }
        [[nodiscard]] size_t getLength() const noexcept { return getRow(); }
    };

    template<class MatrixType>
    class MatrixBlock<MatrixType, Dynamic, Dynamic> : public LValueMatrix<MatrixBlock<MatrixType, Dynamic, Dynamic>> {
    public:
        using Base = LValueMatrix<MatrixBlock<MatrixType, Dynamic, Dynamic>>;
        using typename Base::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        MatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        MatrixBlock(const MatrixBlock&) = default;
        MatrixBlock(MatrixBlock&&) noexcept = delete;
        ~MatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        MatrixBlock& operator=(const MatrixBlock& m) { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        MatrixBlock& operator=(MatrixBlock&& m) noexcept { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        [[nodiscard]] ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t col) const;
        /* Operations */
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] size_t getColumn() const noexcept { return colCount; }
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