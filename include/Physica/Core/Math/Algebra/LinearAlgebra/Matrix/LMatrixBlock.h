/*
 * Copyright 2021-2022 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/ContinuousVector.h"

namespace Physica::Core {
    template<class MatrixType> class RowLVector;
    template<class MatrixType> class ColLVector;
    template<class MatrixType, size_t Row = Dynamic, size_t Column = Dynamic> class LMatrixBlock;

    namespace Internal {
        template<class MatrixType>
        class Traits<RowLVector<MatrixType>> {
            using VectorType = RowLVector<MatrixType>;
            constexpr static bool isRowMatrix = DenseMatrixOption::isRowMatrix<MatrixType>();
        public:
            using Base = typename std::conditional<isRowMatrix, ContinuousVector<VectorType>, LValueVector<VectorType>>::type;
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = Traits<MatrixType>::ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = Traits<MatrixType>::MaxColumnAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };

        template<class MatrixType>
        class Traits<ColLVector<MatrixType>> {
            using VectorType = ColLVector<MatrixType>;
            constexpr static bool isColumnMatrix = DenseMatrixOption::isColumnMatrix<MatrixType>();
        public:
            using Base = typename std::conditional<isColumnMatrix, ContinuousVector<VectorType>, LValueVector<VectorType>>::type;
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = Traits<MatrixType>::RowAtCompile;
            constexpr static size_t MaxSizeAtCompile = Traits<MatrixType>::MaxRowAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };

        template<class MatrixType, size_t Row, size_t Column>
        class Traits<LMatrixBlock<MatrixType, Row, Column>> {
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static int MatrixOption = MatrixType::MatrixOption;
            constexpr static size_t RowAtCompile = Row;
            constexpr static size_t ColumnAtCompile = Column;
            constexpr static size_t MaxRowAtCompile = Row;
            constexpr static size_t MaxColumnAtCompile = Column;
            constexpr static size_t SizeAtCompile = Row * Column;
            constexpr static size_t MaxSizeAtCompile = SizeAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }
    /**
     * \class RowLVector and \class ColLVector is designed to implement \class LMatrixBlock, and they can be used indepently.
     */
    template<class MatrixType>
    class RowLVector : public Internal::Traits<RowLVector<MatrixType>>::Base {
    public:
        using Base = typename Internal::Traits<RowLVector<MatrixType>>::Base;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t row;
        size_t fromCol;
        size_t colCount;
    public:
        RowLVector(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : mat(mat_), row(row_), fromCol(fromCol_), colCount(colCount_) {
            assert(row < mat.getRow());
            assert(fromCol + colCount <= mat.getColumn());
        }
        RowLVector(const RowLVector&) = delete;
        RowLVector(RowLVector&&) noexcept = delete;
        ~RowLVector() = default;
        /* Operators */
        RowLVector& operator=(const RowLVector& v) { v.assignTo(*this); return *this; }
        RowLVector& operator=(RowLVector&& v) noexcept { return operator=(std::cref(v)); }
        using Base::operator=;
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < colCount); return mat(row, fromCol + index); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < colCount); return mat(row, fromCol + index); }
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == colCount); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return colCount; }
    };

    template<class MatrixType>
    class ColLVector : public Internal::Traits<ColLVector<MatrixType>>::Base {
    public:
        using Base = typename Internal::Traits<ColLVector<MatrixType>>::Base;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t col;
    public:
        ColLVector(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : mat(mat_), fromRow(fromRow_), rowCount(rowCount_), col(col_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(col < mat.getColumn());
        }
        ColLVector(const ColLVector&) = delete;
        ColLVector(ColLVector&&) noexcept = delete;
        ~ColLVector() = default;
        /* Operators */
        ColLVector& operator=(const ColLVector& v) { v.assignTo(*this); return *this; }
        ColLVector& operator=(ColLVector&& v) noexcept { return operator=(std::cref(v)); }
        using Base::operator=;
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < rowCount); return mat(fromRow + index, col); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < rowCount); return mat(fromRow + index, col); }
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return rowCount; }
    };

    template<class MatrixType>
    class LMatrixBlock<MatrixType, 1, Dynamic> : public LValueMatrix<LMatrixBlock<MatrixType, 1, Dynamic>>
                                              , public RowLVector<MatrixType> {
    public:
        using Base = LValueMatrix<LMatrixBlock<MatrixType, 1, Dynamic>>;
        using VectorBase = RowLVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        LMatrixBlock(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : VectorBase(mat_, row_, fromCol_, colCount_) {}
        LMatrixBlock(const LMatrixBlock&) = delete;
        LMatrixBlock(LMatrixBlock&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        using VectorBase::operator=;
        [[nodiscard]] ScalarType& operator()([[maybe_unused]] size_t row, size_t col) { assert(row == 0); return VectorBase::operator[](col); }
        [[nodiscard]] const ScalarType& operator()([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return VectorBase::operator[](col); }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == getColumn()); }
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getColumn() const noexcept { return VectorBase::getLength(); }
        using VectorBase::max;
        using VectorBase::min;
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] typename VectorBase::Base& asVector() noexcept { return *this; }
        [[nodiscard]] const typename VectorBase::Base& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class LMatrixBlock<MatrixType, Dynamic, 1> : public LValueMatrix<LMatrixBlock<MatrixType, Dynamic, 1>>
                                              , public ColLVector<MatrixType> {
    public:
        using Base = LValueMatrix<LMatrixBlock<MatrixType, Dynamic, 1>>;
        using VectorBase = ColLVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        LMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        LMatrixBlock(const LMatrixBlock&) = delete;
        LMatrixBlock(LMatrixBlock&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        using VectorBase::operator=;
        [[nodiscard]] ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) { assert(col == 0); return VectorBase::operator[](row); }
        [[nodiscard]] const ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) const { assert(col == 0); return VectorBase::operator[](row); }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow() && col == 1); }
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        [[nodiscard]] size_t getRow() const noexcept { return VectorBase::getLength(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return 1; }
        using VectorBase::max;
        using VectorBase::min;
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] typename VectorBase::Base& asVector() noexcept { return *this; }
        [[nodiscard]] const typename VectorBase::Base& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class LMatrixBlock<MatrixType, Dynamic, Dynamic> : public LValueMatrix<LMatrixBlock<MatrixType, Dynamic, Dynamic>> {
    public:
        using Base = LValueMatrix<LMatrixBlock<MatrixType, Dynamic, Dynamic>>;
        using typename Base::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        LMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        LMatrixBlock(const LMatrixBlock&) = delete;
        LMatrixBlock(LMatrixBlock&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        LMatrixBlock& operator=(const LMatrixBlock& m) { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        LMatrixBlock& operator=(LMatrixBlock&& m) noexcept { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        [[nodiscard]] ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t col) const;
        /* Operations */
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] size_t getColumn() const noexcept { return colCount; }
    };

    template<class MatrixType>
    LMatrixBlock<MatrixType, Dynamic, Dynamic>::LMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_)
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getColumn());
    }

    template<class MatrixType>
    typename LMatrixBlock<MatrixType, Dynamic, Dynamic>::ScalarType&
    LMatrixBlock<MatrixType, Dynamic, Dynamic>::operator()(size_t row, size_t col) {
        assert(row < rowCount);
        assert(col < colCount);
        return mat(row + fromRow, col + fromCol);
    }

    template<class MatrixType>
    const typename LMatrixBlock<MatrixType, Dynamic, Dynamic>::ScalarType&
    LMatrixBlock<MatrixType, Dynamic, Dynamic>::operator()(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat(row + fromRow, col + fromCol);
    }
}