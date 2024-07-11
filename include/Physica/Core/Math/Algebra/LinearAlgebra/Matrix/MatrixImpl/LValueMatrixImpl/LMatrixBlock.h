/*
 * Copyright 2021-2024 WeiBo He.
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

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h>

namespace Physica::Core {
    template<class MatrixType, size_t Row = Dynamic, size_t Column = Dynamic> class LMatrixBlock;
    /**
     * \class RowLVector and \class ColLVector is designed to implement \class LMatrixBlock, and they can be used indepently.
     */
    template<class MatrixType>
    class RowLVector : public LValueVector<RowLVector<MatrixType>> {
    public:
        using Base = LValueVector<RowLVector<MatrixType>>;
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
        RowLVector(const RowLVector&) = default;
        RowLVector(RowLVector&&) noexcept = default;
        ~RowLVector() = default;
        /* Operators */
        RowLVector& operator=(const RowLVector& v) { v.assignTo(*this); return *this; }
        RowLVector& operator=(RowLVector&& v) noexcept { return operator=(std::cref(v)); }
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return colCount; }
        [[nodiscard]] __host__ __device__ inline ScalarType* data_ptr(size_t index) { assert(index < colCount); return mat.data_ptr(row, fromCol + index); }
        [[nodiscard]] __host__ __device__ inline const ScalarType* data_ptr(size_t index) const { assert(index < colCount); return mat.data_ptr(row, fromCol + index); }
    };

    template<class MatrixType>
    class ColLVector : public LValueVector<ColLVector<MatrixType>> {
    public:
        using Base = LValueVector<ColLVector<MatrixType>>;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t col;
        size_t fromRow;
        size_t rowCount;
    public:
        ColLVector(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : mat(mat_), col(col_), fromRow(fromRow_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(col < mat.getColumn());
        }
        ColLVector(const ColLVector&) = default;
        ColLVector(ColLVector&&) noexcept = default;
        ~ColLVector() = default;
        /* Operators */
        ColLVector& operator=(const ColLVector& v) { v.assignTo(*this); return *this; }
        ColLVector& operator=(ColLVector&& v) noexcept { return operator=(std::cref(v)); }
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ inline ScalarType* data_ptr(size_t index) { assert(index < rowCount); return mat.data_ptr(fromRow + index, col); }
        [[nodiscard]] __host__ __device__ inline const ScalarType* data_ptr(size_t index) const { assert(index < rowCount); return mat.data_ptr(fromRow + index, col); }
    };

    template<class MatrixType>
    class LMatrixBlock<MatrixType, 1, Dynamic> : public LValueMatrix<LMatrixBlock<MatrixType, 1, Dynamic>>
                                               , public RowLVector<MatrixType> {
    public:
        using This = LMatrixBlock<MatrixType, 1, Dynamic>;
        using Base = LValueMatrix<This>;
        using VectorBase = RowLVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        LMatrixBlock(MatrixType& mat_, size_t row_, size_t fromCol_, size_t colCount_) : VectorBase(mat_, row_, fromCol_, colCount_) {}
        LMatrixBlock(const LMatrixBlock&) = default;
        LMatrixBlock(LMatrixBlock&&) noexcept = default;
        ~LMatrixBlock() = default;
        /* Operators */
        This& operator=(const ScalarType& s) { VectorBase::operator=(s); return *this; }
        using Base::operator=;
        using VectorBase::operator=;
        [[nodiscard]] ScalarType& operator()([[maybe_unused]] size_t row, size_t col) { assert(row == 0); return VectorBase::operator[](col); }
        [[nodiscard]] const ScalarType& operator()([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return VectorBase::operator[](col); }
        void operator+=(const ScalarType& s) { VectorBase::operator+=(s); }
        void operator-=(const ScalarType& s) { VectorBase::operator-=(s); }
        void operator*=(const ScalarType& s) { VectorBase::operator*=(s); }
        void operator/=(const ScalarType& s) { VectorBase::operator/=(s); }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == getColumn()); }

        using Base::row;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::conjugate;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        [[nodiscard]] constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] size_t getColumn() const noexcept { return VectorBase::getLength(); }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class LMatrixBlock<MatrixType, Dynamic, 1> : public LValueMatrix<LMatrixBlock<MatrixType, Dynamic, 1>>
                                               , public ColLVector<MatrixType> {
    public:
        using This = LMatrixBlock<MatrixType, Dynamic, 1>;
        using Base = LValueMatrix<This>;
        using VectorBase = ColLVector<MatrixType>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        LMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        LMatrixBlock(const LMatrixBlock&) = default;
        LMatrixBlock(LMatrixBlock&&) noexcept = default;
        ~LMatrixBlock() = default;
        /* Operators */
        This& operator=(const ScalarType& s) { VectorBase::operator=(s); return *this; }
        using Base::operator=;
        using VectorBase::operator=;
        void operator+=(const ScalarType& s) { VectorBase::operator+=(s); }
        void operator-=(const ScalarType& s) { VectorBase::operator-=(s); }
        void operator*=(const ScalarType& s) { VectorBase::operator*=(s); }
        void operator/=(const ScalarType& s) { VectorBase::operator/=(s); }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow() && col == 1); }

        using Base::col;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::conjugate;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        [[nodiscard]] size_t getRow() const noexcept { return VectorBase::getLength(); }
        [[nodiscard]] constexpr static size_t getColumn() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t row, [[maybe_unused]] size_t column) { assert(col == 0); return VectorBase::operator[](row); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t row, [[maybe_unused]] size_t column) const { assert(col == 0); return VectorBase::operator[](row); }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class LMatrixBlock<MatrixType, Dynamic, Dynamic> : public LValueMatrix<LMatrixBlock<MatrixType, Dynamic, Dynamic>> {
    public:
        using This = LMatrixBlock<MatrixType, Dynamic, Dynamic>;
        using Base = LValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        LMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        LMatrixBlock(const LMatrixBlock&) = default;
        LMatrixBlock(LMatrixBlock&&) noexcept = default;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        LMatrixBlock& operator=(const LMatrixBlock& m) { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        LMatrixBlock& operator=(LMatrixBlock&& m) noexcept { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        /* Operations */
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] size_t getColumn() const noexcept { return colCount; }
        [[nodiscard]] __host__ __device__ inline ScalarType* data_ptr(size_t row, size_t column);
        [[nodiscard]] __host__ __device__ inline const ScalarType* data_ptr(size_t row, size_t column) const;
        [[nodiscard]] This& asMatrix() noexcept { return *this; }
        [[nodiscard]] const This& asMatrix() const noexcept { return *this; }
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
    __host__ __device__ typename LMatrixBlock<MatrixType, Dynamic, Dynamic>::ScalarType*
    LMatrixBlock<MatrixType, Dynamic, Dynamic>::data_ptr(size_t row, size_t col) {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.data_ptr(row + fromRow, col + fromCol);
    }

    template<class MatrixType>
    __host__ __device__ const typename LMatrixBlock<MatrixType, Dynamic, Dynamic>::ScalarType*
    LMatrixBlock<MatrixType, Dynamic, Dynamic>::data_ptr(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.data_ptr(row + fromRow, col + fromCol);
    }
}

namespace Physica {
    using namespace Core;

    template<class MatrixType>
    class Traits<RowLVector<MatrixType>> {
        constexpr static bool isRowMatrix = MatrixOption::isRowMatrix<MatrixType>();
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = Traits<MatrixType>::ColumnAtCompile;
        constexpr static size_t MaxSizeAtCompile = Traits<MatrixType>::MaxColumnAtCompile;

        constexpr static bool FastAssign = false;
    };

    template<class MatrixType>
    class Traits<ColLVector<MatrixType>> {
        constexpr static bool isColumnMatrix = MatrixOption::isColumnMatrix<MatrixType>();
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = Traits<MatrixType>::RowAtCompile;
        constexpr static size_t MaxSizeAtCompile = Traits<MatrixType>::MaxRowAtCompile;

        constexpr static bool FastAssign = false;
    };

    template<class MatrixType, size_t Row, size_t Column>
    class Traits<LMatrixBlock<MatrixType, Row, Column>> {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static int Option = MatrixType::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColumnAtCompile = Column;
        constexpr static size_t MaxRowAtCompile = Row;
        constexpr static size_t MaxColumnAtCompile = Column;
        constexpr static size_t SizeAtCompile = Row * Column;
        constexpr static size_t MaxSizeAtCompile = SizeAtCompile;
    };
}
