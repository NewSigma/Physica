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

    template<Matrix T, size_t Col>
    class ContinuousMatrixBlock<T, 1, Col> : public ContinuousVector<ContinuousMatrixBlock<T, 1, Col>> {
        using This = ContinuousMatrixBlock<T, 1, Col>;
        using Base = ContinuousVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::isReverseDiff;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        PtrTy pVecHead;
        size_t colCount;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<T>& mat, size_t fromRow, [[maybe_unused]] size_t rowCount_, size_t fromCol, size_t colCount_)
                : ContinuousMatrixBlock(mat, fromRow, fromCol, colCount_) {
            assert(rowCount_ == 1);
        }
        ContinuousMatrixBlock(ContinuousMatrix<T>& mat, size_t fromRow, size_t fromCol, size_t colCount_)
                : pVecHead(mat.data_ptr(fromRow, fromCol)), colCount(colCount_) {
            assert(fromRow < mat.getRow());
            assert(fromCol + colCount <= mat.getCol());
        }
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        [[nodiscard]] This& row(size_t r) {
            assert(r == 0);
            return *this;
        }
        [[nodiscard]] const This& row(size_t r) const {
            assert(r == 0);
            return *this;
        }
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < colCount && "[Error]: Index overflow");
            return pVecHead + index;
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    };

    template<Matrix T, size_t Row>
    class ContinuousMatrixBlock<T, Row, 1> : public ContinuousVector<ContinuousMatrixBlock<T, Row, 1>> {
        using This = ContinuousMatrixBlock<T, Row, 1>;
        using Base = ContinuousVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::isReverseDiff;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        PtrTy pVecHead;
        size_t rowCount;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<T>& mat, size_t fromRow, size_t rowCount_, size_t fromCol_, [[maybe_unused]] size_t colCount_)
                : ContinuousMatrixBlock(mat, fromRow, rowCount_, fromCol_) {
            assert(colCount_ == 1);
        }
        ContinuousMatrixBlock(ContinuousMatrix<T>& mat, size_t fromRow, size_t rowCount_, size_t fromCol)
                : pVecHead(mat.data_ptr(fromRow, fromCol)), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        [[nodiscard]] This& col([[maybe_unused]] size_t c) {
            assert(c == 0);
            return *this;
        }
        [[nodiscard]] const This& col(size_t c) const {
            assert(c == 0);
            return *this;
        }
        void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < rowCount && "[Error]: Index overflow");
            return pVecHead + index;
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    };

    template<Matrix T>
    class ContinuousMatrixBlock<T, 1, 1> : public ContinuousVector<ContinuousMatrixBlock<T, 1, 1>> {
        using This = ContinuousMatrixBlock<T, 1, 1>;
        using Base = ContinuousVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::isReverseDiff;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        PtrTy pVecHead;
        size_t rowCount;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<T>& mat, size_t fromRow, size_t rowCount_, size_t fromCol)
                : pVecHead(mat.data_ptr(fromRow, fromCol)), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        [[nodiscard]] This& row(size_t r) {
            assert(r == 0);
            return *this;
        }
        [[nodiscard]] const This& row(size_t r) const {
            assert(r == 0);
            return *this;
        }
        [[nodiscard]] This& col(size_t c) {
            assert(c == 0);
            return *this;
        }
        [[nodiscard]] const This& col(size_t c) const {
            assert(c == 0);
            return *this;
        }
        void resize([[maybe_unused]] size_t length) { assert(length == 1); }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t index) {
            assert(index == 0 && "[Error]: Index overflow");
            return pVecHead;
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr([[maybe_unused]] size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    };

    template<Matrix T, size_t Row, size_t Col>
    class ContinuousMatrixBlock<T, Row, Col> : public LValueMatrix<ContinuousMatrixBlock<T, Row, Col>> {
        using This = ContinuousMatrixBlock<T, Row, Col>;
        using Base = LValueMatrix<This>;
    public:
        using Base::isComplex;
        using typename Base::ScalarType;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        T& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        ContinuousMatrixBlock(ContinuousMatrix<T>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        ContinuousMatrixBlock(const This&) = default;
        ContinuousMatrixBlock(This&&) noexcept = default;
        ~ContinuousMatrixBlock() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t row, size_t col) const;
    };

    template<Matrix T, size_t Row, size_t Col>
    ContinuousMatrixBlock<T, Row, Col>::ContinuousMatrixBlock(
            ContinuousMatrix<T>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
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

    template<Matrix T, size_t Row, size_t Col>
    __host__ __device__ ContinuousMatrixBlock<T, Row, Col>::PtrTy
    ContinuousMatrixBlock<T, Row, Col>::data_ptr(size_t row, size_t col) {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.data_ptr(row + fromRow, col + fromCol);
    }

    template<Matrix T, size_t Row, size_t Col>
    __host__ __device__ ContinuousMatrixBlock<T, Row, Col>::ConstPtrTy
    ContinuousMatrixBlock<T, Row, Col>::data_ptr(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.data_ptr(row + fromRow, col + fromCol);
    }
}

namespace Physica {
    template<Matrix T, size_t Row, size_t Col>
    class Traits<ContinuousMatrixBlock<T, Row, Col>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = Row * Col;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = true;
    };
}
