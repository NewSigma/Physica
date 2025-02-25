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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.h"

namespace Physica {
    template<class MatrixType, size_t Row = Dynamic, size_t Col = Dynamic> class LMatrixBlock;

    template<Matrix T, size_t Col>
    class LMatrixBlock<T, 1, Col> : public LValueVector<LMatrixBlock<T, 1, Col>> {
        using This = LMatrixBlock<T, 1, Col>;
        using Base = LValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        T& mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        LMatrixBlock(T& mat_, size_t fromRow_, size_t fromCol_, size_t colCount_) : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), colCount(colCount_) {
            assert(fromRow < mat.getRow());
            assert(fromCol + colCount <= mat.getCol());
        }
        LMatrixBlock(const This&) = default;
        LMatrixBlock(This&&) noexcept = default;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < colCount);
            return mat.data_ptr(fromRow, fromCol + index);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    };

    template<Matrix T, size_t Row>
    class LMatrixBlock<T, Row, 1> : public LValueVector<LMatrixBlock<T, Row, 1>> {
        using This = LMatrixBlock<T, Row, 1>;
        using Base = LValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        T& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        LMatrixBlock(T& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        LMatrixBlock(const This&) = default;
        LMatrixBlock(This&&) noexcept = default;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < rowCount);
            return mat.data_ptr(fromRow + index, fromCol);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    };

    template<Matrix T>
    class LMatrixBlock<T, 1, 1> : public LValueVector<LMatrixBlock<T, 1, 1>> {
        using This = LMatrixBlock<T, 1, 1>;
        using Base = LValueVector<This>;
    public:
        using typename Base::ScalarType;
        using typename Base::ValueType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        T& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        LMatrixBlock(T& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(mat_), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        LMatrixBlock(const This&) = delete;
        LMatrixBlock(This&&) noexcept = delete;
        ~LMatrixBlock() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t index) {
            assert(index == 0 && "[Error]: Index overflow");
            return mat.data_ptr(fromRow, fromCol);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr([[maybe_unused]] size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    };

    template<Matrix T>
    class LMatrixBlock<T, Dynamic, Dynamic> : public LValueMatrix<LMatrixBlock<T, Dynamic, Dynamic>> {
        using This = LMatrixBlock<T, Dynamic, Dynamic>;
        using Base = LValueMatrix<This>;
    public:
        using typename Base::ScalarType;
        using Base::isComplex;
        using Base::SizeAtCompile;
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
        LMatrixBlock(T& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
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
    };

    template<Matrix T>
    LMatrixBlock<T, Dynamic, Dynamic>::LMatrixBlock(T& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_)
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<Matrix T>
    __host__ __device__ inline LMatrixBlock<T, Dynamic, Dynamic>::PtrTy
    LMatrixBlock<T, Dynamic, Dynamic>::data_ptr(size_t row, size_t col) {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.data_ptr(row + fromRow, col + fromCol);
    }

    template<Matrix T>
    __host__ __device__ LMatrixBlock<T, Dynamic, Dynamic>::ConstPtrTy
    LMatrixBlock<T, Dynamic, Dynamic>::data_ptr(size_t row, size_t col) const {
        return const_cast<This&>(*this).data_ptr(row, col);
    }
}

namespace Physica {
    template<Matrix T, size_t Row, size_t Col>
    class Traits<LMatrixBlock<T, Row, Col>> {
    public:
        using ScalarType = T::ScalarType;
        constexpr static int Option = T::Option;
        constexpr static size_t RowAtCompile = Row;
        constexpr static size_t ColAtCompile = Col;
        constexpr static size_t SizeAtCompile = Row * Col;

        constexpr static bool FastAssign = false;
        constexpr static bool FastPacket = false;
    };
}
