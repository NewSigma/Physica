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

namespace Physica::Core {
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class ContinuousMatrixBlock;

    template<class MatrixType, size_t Length>
    class RowContinuousVector : public ContinuousVector<RowContinuousVector<MatrixType, Length>> {
        using This = RowContinuousVector<MatrixType, Length>;
        using Base = ContinuousVector<This>;
    public:
        using typename Base::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    private:
        PtrTy pVecHead;
        size_t colCount;
    public:
        RowContinuousVector(ContinuousMatrix<MatrixType>& mat, size_t row, size_t fromCol, size_t colCount_)
                : pVecHead(mat.data_ptr(row, fromCol)), colCount(colCount_) {
            assert(row < mat.getRow());
            assert(fromCol + colCount <= mat.getCol());
        }
        RowContinuousVector(const RowContinuousVector&) = default;
        RowContinuousVector(RowContinuousVector&&) noexcept = default;
        ~RowContinuousVector() = default;
        /* Operators */
        using Base::operator=;
        RowContinuousVector& operator=(const RowContinuousVector& v) { v.assignTo(*this); return *this; }
        RowContinuousVector& operator=(RowContinuousVector&& v) noexcept { return operator=(std::cref(v)); }
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        void makeContinuous() { /* Nothing, it is matrix's responsibility */ }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Length == Dynamic ? colCount : Length; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < colCount && "[Error]: Index overflow");
            return pVecHead + index;
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    private:
        using Base::makeContinuous;
    };

    template<class MatrixType, size_t Length>
    class ColContinuousVector : public ContinuousVector<ColContinuousVector<MatrixType, Length>> {
        using This = ColContinuousVector<MatrixType, Length>;
        using Base = ContinuousVector<This>;
    public:
        using typename Base::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    private:
        PtrTy pVecHead;
        size_t rowCount;
    public:
        ColContinuousVector(ContinuousMatrix<MatrixType>& mat, size_t fromRow, size_t rowCount_, size_t col)
                : pVecHead(mat.data_ptr(fromRow, col)), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(col < mat.getCol());
        }
        ColContinuousVector(const ColContinuousVector&) = default;
        ColContinuousVector(ColContinuousVector&&) noexcept = default;
        ~ColContinuousVector() = default;
        /* Operators */
        using Base::operator=;
        ColContinuousVector& operator=(const ColContinuousVector& v) { v.assignTo(*this); return *this; }
        ColContinuousVector& operator=(ColContinuousVector&& v) noexcept { return operator=(std::cref(v)); }
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        void makeContinuous() { /* Nothing, it is matrix's responsibility */ }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Length == Dynamic ? rowCount : Length; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < rowCount && "[Error]: Index overflow");
            return pVecHead + index;
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    private:
        using Base::makeContinuous;
    };

    template<class MatrixType, size_t Col>
    class ContinuousMatrixBlock<MatrixType, 1, Col>
            : public LValueMatrix<ContinuousMatrixBlock<MatrixType, 1, Col>>
            , public Traits<ContinuousMatrixBlock<MatrixType, 1, Col>>::VectorBase {
        using This = ContinuousMatrixBlock<MatrixType, 1, Col>;
    public:
        using Base = LValueMatrix<This>;
        using Base::isComplex;
        using VectorBase = typename Traits<This>::VectorBase;
        using ScalarType = typename MatrixType::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, [[maybe_unused]] size_t rowCount_, size_t fromCol_, size_t colCount_)
                : ContinuousMatrixBlock(mat_, fromRow_, fromCol_, colCount_) {
            assert(rowCount_ == 1);
        }
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t row_, size_t fromCol_, size_t colCount_)
                : VectorBase(mat_.getDerived(), row_, fromCol_, colCount_) {}
        ContinuousMatrixBlock(const ContinuousMatrixBlock&) = default;
        ContinuousMatrixBlock(ContinuousMatrixBlock&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        ContinuousMatrixBlock& operator=(const ContinuousMatrixBlock& m) { VectorBase::operator=(m.asVector()); return *this; }
        ContinuousMatrixBlock& operator=(ContinuousMatrixBlock&& m) noexcept { VectorBase::operator=(m.asVector()); return *this; }
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
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == getCol()); }

        using VectorBase::format;
        using Base::row;
        using Base::transpose;
        using Base::conjugate;
        using Base::hermite;

        using VectorBase::random_uniform;
        using VectorBase::random_normal;
        using VectorBase::random_any;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        using VectorBase::data_ptr;
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Col == Dynamic ? VectorBase::getLength() : Col; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t row, size_t col) { assert(row == 0); return VectorBase::data_ptr(col); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return VectorBase::data_ptr(col); }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row>
    class ContinuousMatrixBlock<MatrixType, Row, 1>
            : public LValueMatrix<ContinuousMatrixBlock<MatrixType, Row, 1>>
            , public Traits<ContinuousMatrixBlock<MatrixType, Row, 1>>::VectorBase {
        using This = ContinuousMatrixBlock<MatrixType, Row, 1>;
    public:
        using Base = LValueMatrix<This>;
        using Base::isComplex;
        using VectorBase = typename Traits<This>::VectorBase;
        using ScalarType = typename MatrixType::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, [[maybe_unused]] size_t colCount_)
                : ContinuousMatrixBlock(mat_, fromRow_, rowCount_, fromCol_) {
            assert(colCount_ == 1);
        }
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : VectorBase(mat_.getDerived(), fromRow_, rowCount_, col_) {}
        ContinuousMatrixBlock(const ContinuousMatrixBlock&) = default;
        ContinuousMatrixBlock(ContinuousMatrixBlock&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        ContinuousMatrixBlock& operator=(const ContinuousMatrixBlock& m) { VectorBase::operator=(m.asVector()); return *this; }
        ContinuousMatrixBlock& operator=(ContinuousMatrixBlock&& m) noexcept { VectorBase::operator=(m.asVector()); return *this; }
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

        using VectorBase::format;
        using Base::col;
        using Base::transpose;
        using Base::conjugate;
        using Base::hermite;

        using VectorBase::random_uniform;
        using VectorBase::random_normal;
        using VectorBase::random_any;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::conjugate;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        using VectorBase::data_ptr;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Row == Dynamic ? VectorBase::getLength() : Row; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCol() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t row, [[maybe_unused]] size_t col) { assert(col == 0); return VectorBase::data_ptr(row); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t row, [[maybe_unused]] size_t col) const { assert(col == 0); return VectorBase::data_ptr(row); }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class ContinuousMatrixBlock<MatrixType, 1, 1>
            : public LValueMatrix<ContinuousMatrixBlock<MatrixType, 1, 1>>
            , public Traits<ContinuousMatrixBlock<MatrixType, 1, 1>>::VectorBase {
        using This = ContinuousMatrixBlock<MatrixType, 1, 1>;
    public:
        using Base = LValueMatrix<This>;
        using Base::isComplex;
        using VectorBase = typename Traits<This>::VectorBase;
        using ScalarType = typename MatrixType::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : VectorBase(mat_.getDerived(), fromRow_, rowCount_, col_) {}
        ContinuousMatrixBlock(const ContinuousMatrixBlock&) = default;
        ContinuousMatrixBlock(ContinuousMatrixBlock&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        ContinuousMatrixBlock& operator=(const ContinuousMatrixBlock& m) { VectorBase::operator=(m.asVector()); return *this; }
        ContinuousMatrixBlock& operator=(ContinuousMatrixBlock&& m) noexcept { VectorBase::operator=(m.asVector()); return *this; }
        using Base::operator=;
        using VectorBase::operator=;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == 1); }
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::conjugate;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCol() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
            assert(row == 0 && col == 0 && "[Error]: Index overflow");
            return VectorBase::data_ptr(0);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) const {
            assert(row == 0 && col == 0 && "[Error]: Index overflow");
            return VectorBase::data_ptr(0);
        }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row, size_t Col>
    class ContinuousMatrixBlock : public LValueMatrix<ContinuousMatrixBlock<MatrixType, Row, Col>> {
    public:
        using This = ContinuousMatrixBlock<MatrixType, Row, Col>;
        using Base = LValueMatrix<This>;
        using Base::isComplex;
        using typename Base::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        ContinuousMatrixBlock(const This&) = delete;
        ContinuousMatrixBlock(This&&) noexcept = delete;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(static_cast<const RValueMatrix<This>&>(m)); return *this; }
        This& operator=(This&& m) noexcept { Base::operator=(static_cast<RValueMatrix<This>&>(m)); return *this; }
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t row, size_t col) const;
        [[nodiscard]] This& asMatrix() noexcept { return *this; }
        [[nodiscard]] const This& asMatrix() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row, size_t Col>
    ContinuousMatrixBlock<MatrixType, Row, Col>::ContinuousMatrixBlock(
            ContinuousMatrix<MatrixType>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_.getDerived())
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<class MatrixType, size_t Row, size_t Col>
    __host__ __device__ typename ContinuousMatrixBlock<MatrixType, Row, Col>::PtrTy
    ContinuousMatrixBlock<MatrixType, Row, Col>::data_ptr(size_t row, size_t col) {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.data_ptr(row + fromRow, col + fromCol);
    }

    template<class MatrixType, size_t Row, size_t Col>
    __host__ __device__ typename ContinuousMatrixBlock<MatrixType, Row, Col>::ConstPtrTy
    ContinuousMatrixBlock<MatrixType, Row, Col>::data_ptr(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.data_ptr(row + fromRow, col + fromCol);
    }
}

namespace Physica {
    using namespace Core;

    template<class MatrixType, size_t Length>
    class Traits<RowContinuousVector<MatrixType, Length>> {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = Length;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };

    template<class MatrixType, size_t Length>
    class Traits<ColContinuousVector<MatrixType, Length>> {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = Length;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };

    template<class MatrixType, size_t Row, size_t Col>
    class Traits<ContinuousMatrixBlock<MatrixType, Row, Col>> {
        constexpr static bool isElementMatrix = MatrixOption::isElementMatrix<MatrixType>();
        constexpr static bool isRowMatrix = MatrixOption::isRowMatrix<MatrixType>();
        constexpr static bool isColMatrix = MatrixOption::isColumnMatrix<MatrixType>();
        using RowVectorType = typename std::conditional<
                isRowMatrix, RowContinuousVector<MatrixType, Col>, RowLVector<MatrixType>>::type;
        using ColVectorType = typename std::conditional<
                isColMatrix, ColContinuousVector<MatrixType, Row>, ColLVector<MatrixType>>::type;
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static int Option = MatrixType::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColumnAtCompile = Col;
        constexpr static size_t SizeAtCompile = Row * Col;
        using VectorBase = typename std::conditional<Col == 1, ColVectorType, RowVectorType>::type;
    };
}
