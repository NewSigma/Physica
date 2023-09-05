/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    template<class Derived> class ContinuousMatrix;
    template<class MatrixType, size_t Length = Dynamic> class RowContinuousVector;
    template<class MatrixType, size_t Length = Dynamic> class ColContinuousVector;
    template<class MatrixType, size_t Row = Dynamic, size_t Column = Dynamic> class ContinuousMatrixBlock;

    namespace Internal {
        template<class MatrixType, size_t Length>
        class Traits<RowContinuousVector<MatrixType, Length>> {
            using VectorType = RowContinuousVector<MatrixType, Length>;
            constexpr static bool isRowMatrix = MatrixOption::isRowMatrix<MatrixType>();
            constexpr static bool isElementMatrix = MatrixOption::isElementMatrix<MatrixType>();
            constexpr static bool isRowVector = isElementMatrix && MatrixType::RowAtCompile == 1;
        public:
            using Base = typename std::conditional<isRowMatrix || isRowVector, ContinuousVector<VectorType>, LValueVector<VectorType>>::type;
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = Length;
            constexpr static size_t MaxSizeAtCompile = Length;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };

        template<class MatrixType, size_t Length>
        class Traits<ColContinuousVector<MatrixType, Length>> {
            using VectorType = ColContinuousVector<MatrixType, Length>;
            constexpr static bool isColumnMatrix = MatrixOption::isColumnMatrix<MatrixType>();
            constexpr static bool isElementMatrix = MatrixOption::isElementMatrix<MatrixType>();
            constexpr static bool isColumnVector = isElementMatrix && MatrixType::ColumnAtCompile == 1;
        public:
            using Base = typename std::conditional<isColumnMatrix || isColumnVector, ContinuousVector<VectorType>, LValueVector<VectorType>>::type;
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static size_t SizeAtCompile = Length;
            constexpr static size_t MaxSizeAtCompile = Length;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };

        template<class MatrixType, size_t Row, size_t Column>
        class Traits<ContinuousMatrixBlock<MatrixType, Row, Column>> {
            constexpr static bool isElementMatrix = MatrixOption::isElementMatrix<MatrixType>();
        public:
            using ScalarType = typename MatrixType::ScalarType;
            constexpr static int Option = MatrixType::Option;
            constexpr static size_t RowAtCompile = Row;
            constexpr static size_t ColumnAtCompile = Column;
            constexpr static size_t MaxRowAtCompile = Row;
            constexpr static size_t MaxColumnAtCompile = Column;
            constexpr static size_t SizeAtCompile = Row * Column;
            constexpr static size_t MaxSizeAtCompile = SizeAtCompile;
            using PacketType = typename Internal::BestPacket<ScalarType, SizeAtCompile>::Type;
        };
    }

    template<class MatrixType, size_t Length>
    class RowContinuousVector : public Internal::Traits<RowContinuousVector<MatrixType, Length>>::Base {
    public:
        using Base = typename Internal::Traits<RowContinuousVector<MatrixType, Length>>::Base;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t row;
        size_t fromCol;
        size_t colCount;
    public:
        RowContinuousVector(ContinuousMatrix<MatrixType>& mat_, size_t row_, size_t fromCol_, size_t colCount_)
                : mat(mat_.getDerived()), row(row_), fromCol(fromCol_), colCount(colCount_) {
            assert(row < mat.getRow());
            assert(fromCol + colCount <= mat.getColumn());
        }
        RowContinuousVector(const RowContinuousVector&) = delete;
        RowContinuousVector(RowContinuousVector&&) noexcept = delete;
        ~RowContinuousVector() = default;
        /* Operators */
        using Base::operator=;
        RowContinuousVector& operator=(const RowContinuousVector& v) { v.assignTo(*this); return *this; }
        RowContinuousVector& operator=(RowContinuousVector&& v) noexcept { return operator=(std::cref(v)); }
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < getLength()); return mat(row, fromCol + index); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < getLength()); return mat(row, fromCol + index); }
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Length == Dynamic ? colCount : Length; }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t index) { return mat.data_ptr(row, fromCol + index); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t index) const { return mat.data_ptr(row, fromCol + index); }
    };

    template<class MatrixType, size_t Length>
    class ColContinuousVector : public Internal::Traits<ColContinuousVector<MatrixType, Length>>::Base {
    public:
        using Base = typename Internal::Traits<ColContinuousVector<MatrixType, Length>>::Base;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        MatrixType& mat;
        size_t col;
        size_t fromRow;
        size_t rowCount;
    public:
        ColContinuousVector(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : mat(mat_.getDerived()), col(col_), fromRow(fromRow_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(col < mat.getColumn());
        }
        ColContinuousVector(const ColContinuousVector&) = delete;
        ColContinuousVector(ColContinuousVector&&) noexcept = delete;
        ~ColContinuousVector() = default;
        /* Operators */
        using Base::operator=;
        ColContinuousVector& operator=(const ColContinuousVector& v) { v.assignTo(*this); return *this; }
        ColContinuousVector& operator=(ColContinuousVector&& v) noexcept { return operator=(std::cref(v)); }
        [[nodiscard]] ScalarType& operator[](size_t index) { assert(index < getLength()); return mat(fromRow + index, col); }
        [[nodiscard]] const ScalarType& operator[](size_t index) const { assert(index < getLength()); return mat(fromRow + index, col); }
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Length == Dynamic ? rowCount : Length; }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t index) { return mat.data_ptr(fromRow + index, col); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t index) const { return mat.data_ptr(fromRow + index, col); }
    };

    template<class MatrixType, size_t Column>
    class ContinuousMatrixBlock<MatrixType, 1, Column> : public LValueMatrix<ContinuousMatrixBlock<MatrixType, 1, Column>>
                                                       , public RowContinuousVector<MatrixType, Column> {
        using This = ContinuousMatrixBlock<MatrixType, 1, Column>;
    public:
        using Base = LValueMatrix<This>;
        using VectorBase = RowContinuousVector<MatrixType, Column>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, [[maybe_unused]] size_t rowCount_, size_t fromCol_, size_t colCount_)
                : ContinuousMatrixBlock(mat_, fromRow_, fromCol_, colCount_) {
            assert(rowCount_ == 1);
        }
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t row_, size_t fromCol_, size_t colCount_) : VectorBase(mat_, row_, fromCol_, colCount_) {}
        ContinuousMatrixBlock(const ContinuousMatrixBlock&) = delete;
        ContinuousMatrixBlock(ContinuousMatrixBlock&&) noexcept = delete;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        ContinuousMatrixBlock& operator=(const ContinuousMatrixBlock& m) { VectorBase::operator=(m.asVector()); return *this; }
        ContinuousMatrixBlock& operator=(ContinuousMatrixBlock&& m) noexcept { VectorBase::operator=(m.asVector()); return *this; }
        template<class T> This& operator=(const ScalarBase<T>& s) { VectorBase::operator=(s); return *this; }
        using Base::operator=;
        using VectorBase::operator=;
        [[nodiscard]] ScalarType& operator()([[maybe_unused]] size_t row, size_t col) { assert(row == 0); return VectorBase::operator[](col); }
        [[nodiscard]] const ScalarType& operator()([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return VectorBase::operator[](col); }
        template<class T> void operator+=(const ScalarBase<T>& s) { VectorBase::operator+=(s); }
        template<class T> void operator-=(const ScalarBase<T>& s) { VectorBase::operator-=(s); }
        template<class T> void operator*=(const ScalarBase<T>& s) { VectorBase::operator*=(s); }
        template<class T> void operator/=(const ScalarBase<T>& s) { VectorBase::operator/=(s); }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == getColumn()); }

        using Base::row;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return Column == Dynamic ? VectorBase::getLength() : Column; }
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] typename VectorBase::Base& asVector() noexcept { return *this; }
        [[nodiscard]] const typename VectorBase::Base& asVector() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row>
    class ContinuousMatrixBlock<MatrixType, Row, 1> : public LValueMatrix<ContinuousMatrixBlock<MatrixType, Row, 1>>
                                                    , public ColContinuousVector<MatrixType, Row> {
        using This = ContinuousMatrixBlock<MatrixType, Row, 1>;
    public:
        using Base = LValueMatrix<This>;
        using VectorBase = ColContinuousVector<MatrixType, Row>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, [[maybe_unused]] size_t colCount_)
                : ContinuousMatrixBlock(mat_, fromRow_, rowCount_, fromCol_) {
            assert(colCount_ == 1);
        }
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t col_) : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        ContinuousMatrixBlock(const ContinuousMatrixBlock&) = delete;
        ContinuousMatrixBlock(ContinuousMatrixBlock&&) noexcept = delete;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        ContinuousMatrixBlock& operator=(const ContinuousMatrixBlock& m) { VectorBase::operator=(m.asVector()); return *this; }
        ContinuousMatrixBlock& operator=(ContinuousMatrixBlock&& m) noexcept { VectorBase::operator=(m.asVector()); return *this; }
        template<class T> This& operator=(const ScalarBase<T>& s) { VectorBase::operator=(s); return *this; }
        using Base::operator=;
        using VectorBase::operator=;
        [[nodiscard]] ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) { assert(col == 0); return VectorBase::operator[](row); }
        [[nodiscard]] const ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) const { assert(col == 0); return VectorBase::operator[](row); }
        template<class T> void operator+=(const ScalarBase<T>& s) { VectorBase::operator+=(s); }
        template<class T> void operator-=(const ScalarBase<T>& s) { VectorBase::operator-=(s); }
        template<class T> void operator*=(const ScalarBase<T>& s) { VectorBase::operator*=(s); }
        template<class T> void operator/=(const ScalarBase<T>& s) { VectorBase::operator/=(s); }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow() && col == 1); }

        using Base::col;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Row == Dynamic ? VectorBase::getLength() : Row; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getColumn() noexcept { return 1; }
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] typename VectorBase::Base& asVector() noexcept { return *this; }
        [[nodiscard]] const typename VectorBase::Base& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class ContinuousMatrixBlock<MatrixType, 1, 1> : public LValueMatrix<ContinuousMatrixBlock<MatrixType, 1, 1>>
                                                  , public ColContinuousVector<MatrixType, 1> {
    public:
        using Base = LValueMatrix<ContinuousMatrixBlock<MatrixType, 1, 1>>;
        using VectorBase = ColContinuousVector<MatrixType, 1>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        ContinuousMatrixBlock(const ContinuousMatrixBlock&) = delete;
        ContinuousMatrixBlock(ContinuousMatrixBlock&&) noexcept = delete;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        ContinuousMatrixBlock& operator=(const ContinuousMatrixBlock& m) { VectorBase::operator=(m.asVector()); return *this; }
        ContinuousMatrixBlock& operator=(ContinuousMatrixBlock&& m) noexcept { VectorBase::operator=(m.asVector()); return *this; }
        using Base::operator=;
        using VectorBase::operator=;
        [[nodiscard]] ScalarType& operator()([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
            assert(row == 0 && col == 0);
            return VectorBase::operator[](0);
        }
        [[nodiscard]] const ScalarType& operator()([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) const {
            assert(row == 0 && col == 0);
            return VectorBase::operator[](0);
        }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == 1); }
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getColumn() noexcept { return 1; }
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

    template<class MatrixType, size_t Row, size_t Column>
    class ContinuousMatrixBlock : public LValueMatrix<ContinuousMatrixBlock<MatrixType, Row, Column>> {
    public:
        using This = ContinuousMatrixBlock<MatrixType, Row, Column>;
        using Base = LValueMatrix<This>;
        using typename Base::ScalarType;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        ContinuousMatrixBlock(const ContinuousMatrixBlock&) = delete;
        ContinuousMatrixBlock(ContinuousMatrixBlock&&) noexcept = delete;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        ContinuousMatrixBlock& operator=(const ContinuousMatrixBlock& m) { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        ContinuousMatrixBlock& operator=(ContinuousMatrixBlock&& m) noexcept { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        [[nodiscard]] ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] const ScalarType& operator()(size_t row, size_t col) const;
        /* Operations */
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return Column == Dynamic ? colCount : Column; }
        [[nodiscard]] This& asMatrix() noexcept { return *this; }
        [[nodiscard]] const This& asMatrix() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row, size_t Column>
    ContinuousMatrixBlock<MatrixType, Row, Column>::ContinuousMatrixBlock(
            ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_.getDerived())
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getColumn());
    }

    template<class MatrixType, size_t Row, size_t Column>
    typename ContinuousMatrixBlock<MatrixType, Row, Column>::ScalarType&
    ContinuousMatrixBlock<MatrixType, Row, Column>::operator()(size_t row, size_t col) {
        assert(row < rowCount);
        assert(col < colCount);
        return mat(row + fromRow, col + fromCol);
    }

    template<class MatrixType, size_t Row, size_t Column>
    const typename ContinuousMatrixBlock<MatrixType, Row, Column>::ScalarType&
    ContinuousMatrixBlock<MatrixType, Row, Column>::operator()(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat(row + fromRow, col + fromCol);
    }
}
