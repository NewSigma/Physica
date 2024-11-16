/*
 * Copyright 2021-2024 Weibo He.
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
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class LMatrixBlock;

    template<class MatrixType, size_t Col>
    class LMatrixBlock<MatrixType, 1, Col> : public LValueMatrix<LMatrixBlock<MatrixType, 1, Col>>
                                           , public LValueVector<LMatrixBlock<MatrixType, 1, Col>> {
        using This = LMatrixBlock<MatrixType, 1, Col>;
        using Base = LValueMatrix<This>;
        using VectorBase = LValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        LMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t fromCol_, size_t colCount_) : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), colCount(colCount_) {
            assert(fromRow < mat.getRow());
            assert(fromCol + colCount <= mat.getCol());
        }
        LMatrixBlock(const This&) = delete;
        LMatrixBlock(This&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        using VectorBase::operator=;
        using VectorBase::operator+=;
        using VectorBase::operator-=;
        using VectorBase::operator*=;
        using VectorBase::operator/=;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t length) { assert(length == colCount); }
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == getCol()); }

        using VectorBase::format;

        using Base::row;
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
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t row, size_t col) {
            assert(row == 0);
            return data_ptr(col);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t row, size_t col) const {
            return const_cast<This&>(*this).data_ptr(row, col);
        }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < colCount);
            return mat.data_ptr(fromRow, fromCol + index);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row>
    class LMatrixBlock<MatrixType, Row, 1> : public LValueMatrix<LMatrixBlock<MatrixType, Row, 1>>
                                           , public LValueVector<LMatrixBlock<MatrixType, Row, 1>> {
        using This = LMatrixBlock<MatrixType, Row, 1>;
        using Base = LValueMatrix<This>;
        using VectorBase = LValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        LMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        LMatrixBlock(const This&) = delete;
        LMatrixBlock(This&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        using VectorBase::operator=;
        using VectorBase::operator+=;
        using VectorBase::operator-=;
        using VectorBase::operator*=;
        using VectorBase::operator/=;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow() && col == 1); }

        using VectorBase::format;

        using Base::col;
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
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCol() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t row, [[maybe_unused]] size_t col) {
            assert(col == 0);
            return data_ptr(row);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t row, size_t col) const {
            return const_cast<This&>(*this).data_ptr(row, col);
        }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < rowCount);
            return mat.data_ptr(fromRow + index, fromCol);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class LMatrixBlock<MatrixType, 1, 1> : public LValueMatrix<LMatrixBlock<MatrixType, 1, 1>>
                                         , public LValueVector<LMatrixBlock<MatrixType, 1, 1>> {
        using This = LMatrixBlock<MatrixType, 1, 1>;
        using Base = LValueMatrix<This>;
        using VectorBase = LValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    private:
        MatrixType& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        LMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        LMatrixBlock(const This&) = delete;
        LMatrixBlock(This&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        using VectorBase::operator=;
        using VectorBase::operator+=;
        using VectorBase::operator-=;
        using VectorBase::operator*=;
        using VectorBase::operator/=;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow() && col == 1); }

        using VectorBase::format;

        using Base::col;
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
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCol() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
            assert(row == 0 && col == 0);
            return data_ptr(0);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t row, size_t col) const {
            return const_cast<This&>(*this).data_ptr(row, col);
        }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t index) {
            assert(index == 0 && "[Error]: Index overflow");
            return mat.data_ptr(fromRow, fromCol);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr([[maybe_unused]] size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
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
        using This = LMatrixBlock<MatrixType, Dynamic, Dynamic>;
        using Base = LValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using Base::isComplex;
        using Base::SizeAtCompile;
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
        LMatrixBlock(MatrixType& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        LMatrixBlock(const This&) = delete;
        LMatrixBlock(This&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return colCount; }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t row, size_t col) const;
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
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<class MatrixType>
    __host__ __device__ inline typename LMatrixBlock<MatrixType, Dynamic, Dynamic>::PtrTy
    LMatrixBlock<MatrixType, Dynamic, Dynamic>::data_ptr(size_t row, size_t col) {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.data_ptr(row + fromRow, col + fromCol);
    }

    template<class MatrixType>
    __host__ __device__ typename LMatrixBlock<MatrixType, Dynamic, Dynamic>::ConstPtrTy
    LMatrixBlock<MatrixType, Dynamic, Dynamic>::data_ptr(size_t row, size_t col) const {
        return const_cast<This&>(*this).data_ptr(row, col);
    }
}

namespace Physica {
    template<class MatrixType, size_t Row, size_t Col>
    class Traits<Core::LMatrixBlock<MatrixType, Row, Col>> {
    public:
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static int Option = MatrixType::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = Row * Col;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
