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

#include "ContinuousMatrixBlock.h"

namespace Physica::Core {
    template<class MatrixType, size_t Col>
    class device_obj<ContinuousMatrixBlock<MatrixType, 1, Col>>
            : public device_obj<LValueMatrix<ContinuousMatrixBlock<MatrixType, 1, Col>>>
            , public device_obj<ContinuousVector<ContinuousMatrixBlock<MatrixType, 1, Col>>> {
        using host_obj = ContinuousMatrixBlock<MatrixType, 1, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
        using VectorBase = device_obj<ContinuousVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using Base::isComplex;
        using Base::isReverseDiff;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    private:
        device_obj<MatrixType>& mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat, size_t fromRow_, [[maybe_unused]] size_t rowCount, size_t fromCol_, size_t colCount_)
                : device_obj(mat, fromRow_, fromCol_, colCount_) {
            assert(rowCount == 1);
        }
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat, size_t fromRow_, size_t fromCol_, size_t colCount_)
                : mat(mat.getDerived()), fromRow(fromRow_), fromCol(fromCol_), colCount(colCount_) {
            assert(fromRow < mat.getRow());
            assert(fromCol + colCount <= mat.getCol());
        }
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
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
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == getCol()); }

        using Base::row;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        using VectorBase::data_ptr;
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return getLength(); }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t row, size_t col) { assert(row == 0); return VectorBase::data_ptr(col); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return VectorBase::data_ptr(col); }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) { return mat.data_ptr(fromRow, fromCol + index); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const { return mat.data_ptr(fromRow, fromCol + index); }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] __host__ __device__ Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row>
    class device_obj<ContinuousMatrixBlock<MatrixType, Row, 1>>
            : public device_obj<LValueMatrix<ContinuousMatrixBlock<MatrixType, Row, 1>>>
            , public device_obj<ContinuousVector<ContinuousMatrixBlock<MatrixType, Row, 1>>> {
        using host_obj = ContinuousMatrixBlock<MatrixType, Row, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
        using VectorBase = device_obj<ContinuousVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using Base::isComplex;
        using Base::isReverseDiff;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    private:
        device_obj<MatrixType>& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat, size_t fromRow_, size_t rowCount_, size_t fromCol_, [[maybe_unused]] size_t colCount)
                : device_obj(mat, fromRow_, rowCount_, fromCol_) {
            assert(colCount == 1);
        }
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(mat.getDerived()), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
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
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow() && col == 1); }

        using Base::col;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        using VectorBase::data_ptr;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return getLength(); }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCol() noexcept { return 1; }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t row, [[maybe_unused]] size_t col) { assert(col == 0); return VectorBase::data_ptr(row); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t row, [[maybe_unused]] size_t col) const { assert(col == 0); return VectorBase::data_ptr(row); }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) { return mat.data_ptr(fromRow + index, fromCol); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const { return mat.data_ptr(fromRow + index, fromCol); }
        /**
         * There are some common functions shared by vector and matrix, it is necessary to decide which function to call explicitly.
         */
        [[nodiscard]] __host__ __device__ Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType>
    class device_obj<ContinuousMatrixBlock<MatrixType, 1, 1>>
            : public device_obj<LValueMatrix<ContinuousMatrixBlock<MatrixType, 1, 1>>>
            , public device_obj<ContinuousVector<ContinuousMatrixBlock<MatrixType, 1, 1>>> {
        using host_obj = ContinuousMatrixBlock<MatrixType, 1, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
        using VectorBase = device_obj<ContinuousVector<host_obj>>;
    public:
        using typename Base::ScalarType;
    protected:
        using PtrTy = typename ScalarType::PtrTy;
        using ConstPtrTy = typename ScalarType::ConstPtrTy;
    private:
        device_obj<MatrixType>& mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(mat.getDerived()), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        using VectorBase::operator=;
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == 1); }
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::max;
        using VectorBase::min;
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getCol() noexcept { return 1; }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
            assert(row == 0 && col == 0 && "[Error]: Index overflow");
            return data_ptr(0);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) const {
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
        [[nodiscard]] __host__ __device__ Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row, size_t Col>
    class device_obj<ContinuousMatrixBlock<MatrixType, Row, Col>>
            : public device_obj<LValueMatrix<ContinuousMatrixBlock<MatrixType, Row, Col>>> {
        using host_obj = ContinuousMatrixBlock<MatrixType, Row, Col>;
        using This = device_obj<host_obj>;
    public:
        using Base = device_obj<LValueMatrix<host_obj>>;
        using typename Base::ScalarType;
    private:
        device_obj<MatrixType>& mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ This& asMatrix() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const This& asMatrix() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row, size_t Col>
    __host__ __device__ device_obj<ContinuousMatrixBlock<MatrixType, Row, Col>>::device_obj(
            device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_.getDerived())
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }
}

namespace Physica {
    template<class MatrixType, size_t Row, size_t Col>
    class Traits<Core::device_obj<Core::ContinuousMatrixBlock<MatrixType, Row, Col>>> : public Traits<Core::ContinuousMatrixBlock<MatrixType, Row, Col>> {};
}
