/*
 * Copyright 2022 WeiBo He.
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
    template<class Derived> class RValueMatrix;
    template<class MatrixType> class RowRVector;
    template<class MatrixType> class ColRVector;
    template<class MatrixType, size_t Row = Dynamic, size_t Column = Dynamic> class RMatrixBlock;

    namespace Internal {
        template<class MatrixType>
        class Traits<RowRVector<MatrixType>> {
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = Traits<MatrixType>::ColumnAtCompile;
            constexpr static size_t MaxSizeAtCompile = Traits<MatrixType>::MaxColumnAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };

        template<class MatrixType>
        class Traits<ColRVector<MatrixType>> {
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = Traits<MatrixType>::RowAtCompile;
            constexpr static size_t MaxSizeAtCompile = Traits<MatrixType>::MaxRowAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };

        template<class MatrixType, size_t Row, size_t Column>
        class Traits<RMatrixBlock<MatrixType, Row, Column>> {
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
     * \class RowRVector and \class ColRVector is designed to implement \class RMatrixBlock, and they can be used indepently.
     */
    template<class MatrixType>
    class RowRVector : public RValueVector<RowRVector<MatrixType>> {
    public:
        using Base = RValueVector<RowRVector<MatrixType>>;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t row;
        size_t fromCol;
        size_t colCount;
    public:
        RowRVector(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : mat(mat_), row(row_), fromCol(fromCol_), colCount(colCount_) {
            assert(row < mat.getRow());
            assert(fromCol + colCount <= mat.getColumn());
        }
        RowRVector(const RowRVector&) = delete;
        RowRVector(RowRVector&&) noexcept = delete;
        ~RowRVector() = default;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { assert(index < colCount); return mat.calc(row, fromCol + index); }
        [[nodiscard]] size_t getLength() const noexcept { return colCount; }
    };

    template<class MatrixType>
    class ColRVector : public RValueVector<ColRVector<MatrixType>> {
    public:
        using Base = RValueVector<ColRVector<MatrixType>>;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t col;
    public:
        ColRVector(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : mat(mat_), fromRow(fromRow_), rowCount(rowCount_), col(col_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(col < mat.getColumn());
        }
        ColRVector(const ColRVector&) = delete;
        ColRVector(ColRVector&&) noexcept = delete;
        ~ColRVector() = default;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { assert(index < rowCount); return mat.calc(fromRow + index, col); }
        [[nodiscard]] size_t getLength() const noexcept { return rowCount; }
    };

    template<class MatrixType>
    class RMatrixBlock<MatrixType, 1, Dynamic> : public RValueMatrix<RMatrixBlock<MatrixType, 1, Dynamic>>
                                               , public RowRVector<MatrixType> {
    public:
        using Base = RValueMatrix<RMatrixBlock<MatrixType, 1, Dynamic>>;
        using VectorBase = RowRVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        RMatrixBlock(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : VectorBase(mat_, row_, fromCol_, colCount_) {}
        RMatrixBlock(const RMatrixBlock&) = delete;
        RMatrixBlock(RMatrixBlock&&) noexcept = delete;
        ~RMatrixBlock() = default;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        /* Getters */
        [[nodiscard]] ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return VectorBase::calc(col); }
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
    class RMatrixBlock<MatrixType, Dynamic, 1> : public RValueMatrix<RMatrixBlock<MatrixType, Dynamic, 1>>
                                              , public ColRVector<MatrixType> {
    public:
        using Base = RValueMatrix<RMatrixBlock<MatrixType, Dynamic, 1>>;
        using VectorBase = ColRVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        RMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        RMatrixBlock(const RMatrixBlock&) = delete;
        RMatrixBlock(RMatrixBlock&&) noexcept = delete;
        ~RMatrixBlock() = default;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, [[maybe_unused]] size_t col) const { assert(col == 0); return VectorBase::calc(row); }
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
    class RMatrixBlock<MatrixType, Dynamic, Dynamic> : public RValueMatrix<RMatrixBlock<MatrixType, Dynamic, Dynamic>> {
    public:
        using Base = RValueMatrix<RMatrixBlock<MatrixType, Dynamic, Dynamic>>;
        using typename Base::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        RMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        RMatrixBlock(const RMatrixBlock&) = delete;
        RMatrixBlock(RMatrixBlock&&) noexcept = delete;
        ~RMatrixBlock() = default;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] size_t getColumn() const noexcept { return colCount; }
    };

    template<class MatrixType>
    RMatrixBlock<MatrixType, Dynamic, Dynamic>::RMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_)
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getColumn());
    }

    template<class MatrixType>
    typename RMatrixBlock<MatrixType, Dynamic, Dynamic>::ScalarType
    RMatrixBlock<MatrixType, Dynamic, Dynamic>::calc(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.calc(row + fromRow, col + fromCol);
    }
}
