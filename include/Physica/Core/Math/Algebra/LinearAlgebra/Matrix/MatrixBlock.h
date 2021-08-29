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
    template<class MatrixType> class RowVector;

    template<class MatrixType> class ColVector;

    template<class MatrixType, size_t Row = Dynamic, size_t Column = Dynamic>
    class MatrixBlock;

    namespace Internal {
        template<class MatrixType>
        class Traits<RowVector<MatrixType>> {
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = Traits<MatrixType>::ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = Traits<MatrixType>::MaxColumnAtCompile;
        };

        template<class MatrixType>
        class Traits<ColVector<MatrixType>> {
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = Traits<MatrixType>::RowAtCompile;
            constexpr static size_t MaxSizeAtCompile = Traits<MatrixType>::MaxRowAtCompile;
        };

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
    /**
     * \class RowVector and \class ColVector is designed to implement \class MatrixBlock, and they can be used indepently.
     */
    template<class MatrixType>
    class RowVector : public LValueVector<RowVector<MatrixType>> {
    public:
        using Base = LValueVector<RowVector<MatrixType>>;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t row;
        size_t fromCol;
        size_t colCount;
    public:
        RowVector(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : mat(mat_), row(row_), fromCol(fromCol_), colCount(colCount_) {
            assert(row < mat.getRow());
            assert(fromCol + colCount <= mat.getRow());
        }
        RowVector(const RowVector&) = default;
        RowVector(RowVector&&) noexcept = delete;
        ~RowVector() = default;
        /* Operators */
        using Base::operator=;
        RowVector& operator=(const RowVector& m) { Base::operator=(m); return *this; }
        RowVector& operator=(RowVector&& m) noexcept { Base::operator=(m); return *this; }
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < colCount); return mat(row, fromCol + index); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < colCount); return mat(row, fromCol + index); }
        /* Operations */
        void resize(size_t length) { assert(length == colCount); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return colCount; }
    };

    template<class MatrixType>
    class ColVector : public LValueVector<ColVector<MatrixType>> {
    public:
        using Base = LValueVector<ColVector<MatrixType>>;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t col;
    public:
        ColVector(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : mat(mat_), fromRow(fromRow_), rowCount(rowCount_), col(col_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(col < mat.getColumn());
        }
        ColVector(const ColVector&) = default;
        ColVector(ColVector&&) noexcept = delete;
        ~ColVector() = default;
        /* Operators */
        using Base::operator=;
        ColVector& operator=(const ColVector& m) { Base::operator=(m); return *this; }
        ColVector& operator=(ColVector&& m) noexcept { Base::operator=(m); return *this; }
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < rowCount); return mat(fromRow + index, col); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < rowCount); return mat(fromRow + index, col); }
        /* Operations */
        void resize(size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return rowCount; }
    };

    template<class MatrixType>
    class MatrixBlock<MatrixType, 1, Dynamic> : public LValueMatrix<MatrixBlock<MatrixType, 1, Dynamic>>
                                              , public RowVector<MatrixType> {
    public:
        using Base = LValueMatrix<MatrixBlock<MatrixType, 1, Dynamic>>;
        using VectorBase = RowVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        MatrixBlock(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : VectorBase(mat_, row_, fromCol_, colCount_) {}
        MatrixBlock(const MatrixBlock&) = default;
        MatrixBlock(MatrixBlock&&) noexcept = delete;
        ~MatrixBlock() = default;
        /* Operators */
        using VectorBase::operator=;
        MatrixBlock& operator=(const MatrixBlock& m) { VectorBase::operator=(m); return *this; }
        MatrixBlock& operator=(MatrixBlock&& m) noexcept { VectorBase::operator=(m); return *this; }
        [[nodiscard]] ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) { assert(row == 0); return VectorBase::operator[](col); }
        [[nodiscard]] const ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) const { assert(row == 0); return VectorBase::operator[](col); }
        /* Operations */
        void resize(size_t row, size_t col) { assert(row == 1 && col == getColumn()); }
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getColumn() const noexcept { return VectorBase::getLength(); }
    };

    template<class MatrixType>
    class MatrixBlock<MatrixType, Dynamic, 1> : public LValueMatrix<MatrixBlock<MatrixType, Dynamic, 1>>
                                              , public ColVector<MatrixType> {
    public:
        using Base = LValueMatrix<MatrixBlock<MatrixType, Dynamic, 1>>;
        using VectorBase = ColVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        MatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        MatrixBlock(const MatrixBlock&) = default;
        MatrixBlock(MatrixBlock&&) noexcept = delete;
        ~MatrixBlock() = default;
        /* Operators */
        using VectorBase::operator=;
        MatrixBlock& operator=(const MatrixBlock& m) { VectorBase::operator=(m); return *this; }
        MatrixBlock& operator=(MatrixBlock&& m) noexcept { VectorBase::operator=(m); return *this; }
        [[nodiscard]] ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) { assert(col == 0); return VectorBase::operator[](row); }
        [[nodiscard]] const ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) const { assert(col == 0); return VectorBase::operator[](row); }
        /* Operations */
        void resize(size_t row, size_t col) { assert(row == getRow() && col == 1); }
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        [[nodiscard]] size_t getRow() const noexcept { return VectorBase::getLength(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return 1; }
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