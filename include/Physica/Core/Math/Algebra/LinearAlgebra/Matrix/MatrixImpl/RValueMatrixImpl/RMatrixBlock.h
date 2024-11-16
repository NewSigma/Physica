/*
 * Copyright 2022-2024 Weibo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/ContinuousVector.h>

namespace Physica::Core {
    template<class Derived> class RValueMatrix;
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class RMatrixBlock;
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
            assert(fromCol + colCount <= mat.getCol());
        }
        RowRVector(const RowRVector&) = delete;
        RowRVector(RowRVector&&) noexcept = delete;
        ~RowRVector() = default;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { assert(index < colCount); return mat.calc(row, fromCol + index); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return colCount; }
    };

    template<class MatrixType>
    class ColRVector : public RValueVector<ColRVector<MatrixType>> {
    public:
        using Base = RValueVector<ColRVector<MatrixType>>;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t col;
        size_t fromRow;
        size_t rowCount;
    public:
        ColRVector(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : mat(mat_), col(col_), fromRow(fromRow_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(col < mat.getCol());
        }
        ColRVector(const ColRVector&) = delete;
        ColRVector(ColRVector&&) noexcept = delete;
        ~ColRVector() = default;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t index) const { assert(index < rowCount); return mat.calc(fromRow + index, col); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return rowCount; }
    };

    template<class MatrixType>
    class RMatrixBlock<MatrixType, 1, Dynamic> : public RValueMatrix<RMatrixBlock<MatrixType, 1, Dynamic>>
                                               , public RowRVector<MatrixType> {
        using This = RMatrixBlock<MatrixType, 1, Dynamic>;
        using Base = RValueMatrix<This>;
    public:
        using VectorBase = RowRVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        RMatrixBlock(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : VectorBase(mat_, row_, fromCol_, colCount_) {}
        RMatrixBlock(const This&) = delete;
        RMatrixBlock(This&&) noexcept = delete;
        ~RMatrixBlock() = default;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;

        using VectorBase::format;
        /* Getters */
        using VectorBase::calc;
        using VectorBase::conjugate;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        [[nodiscard]] ScalarType calc([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return VectorBase::calc(col); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return VectorBase::getLength(); }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class RMatrixBlock<MatrixType, Dynamic, 1> : public RValueMatrix<RMatrixBlock<MatrixType, Dynamic, 1>>
                                               , public ColRVector<MatrixType> {
        using This = RMatrixBlock<MatrixType, Dynamic, 1>;
        using Base = RValueMatrix<This>;
    public:
        using VectorBase = ColRVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        RMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        RMatrixBlock(const This&) = delete;
        RMatrixBlock(This&&) noexcept = delete;
        ~RMatrixBlock() = default;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;

        using VectorBase::format;
        /* Getters */
        using VectorBase::calc;
        using VectorBase::conjugate;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        [[nodiscard]] ScalarType calc(size_t row, [[maybe_unused]] size_t col) const { assert(col == 0); return VectorBase::calc(row); }
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return VectorBase::getLength(); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCol() noexcept { return 1; }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class RMatrixBlock<MatrixType, Dynamic, Dynamic> : public RValueMatrix<RMatrixBlock<MatrixType, Dynamic, Dynamic>> {
        using This = RMatrixBlock<MatrixType, Dynamic, Dynamic>;
        using Base = RValueMatrix<This>;
    public:
        using typename Base::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        RMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        RMatrixBlock(const This&) = delete;
        RMatrixBlock(This&&) noexcept = delete;
        ~RMatrixBlock() = default;
        /* Getters */
        [[nodiscard]] ScalarType calc(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return colCount; }
    };

    template<class MatrixType>
    RMatrixBlock<MatrixType, Dynamic, Dynamic>::RMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_)
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<class MatrixType>
    typename RMatrixBlock<MatrixType, Dynamic, Dynamic>::ScalarType
    RMatrixBlock<MatrixType, Dynamic, Dynamic>::calc(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.calc(row + fromRow, col + fromCol);
    }
}

namespace Physica {
    using namespace Core;

    template<class MatrixType>
    class Traits<RowRVector<MatrixType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = Traits<MatrixType>::ColAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };

    template<class MatrixType>
    class Traits<ColRVector<MatrixType>> {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = Traits<MatrixType>::RowAtCompile;
        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };

    template<class MatrixType, size_t Row, size_t Col>
    class Traits<RMatrixBlock<MatrixType, Row, Col>> {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static int Option = MatrixType::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = Row * Col;
    };
}
