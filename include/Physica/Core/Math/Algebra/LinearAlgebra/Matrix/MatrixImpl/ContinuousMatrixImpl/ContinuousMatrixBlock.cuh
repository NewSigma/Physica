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
    template<class MatrixType, size_t Length>
    class device_obj<RowContinuousVector<MatrixType, Length>> : public Traits<device_obj<RowContinuousVector<MatrixType, Length>>>::Base {
    public:
        using host_obj = RowContinuousVector<MatrixType, Length>;
        using This = device_obj<host_obj>;
        using Base = typename Traits<This>::Base;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        device_obj<MatrixType>& mat;
        size_t row;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t row_, size_t fromCol_, size_t colCount_)
                : mat(mat_.getDerived()), row(row_), fromCol(fromCol_), colCount(colCount_) {
            assert(row < mat.getRow());
            assert(fromCol + colCount <= mat.getColumn());
        }
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        __device__ device_obj& operator=(const device_obj& v) { v.assignTo(*this); return *this; }
        __device__ device_obj& operator=(device_obj&& v) noexcept { return operator=(std::cref(v)); }
        [[nodiscard]] __device__ ScalarType& operator[](size_t index) { assert(index < getLength()); return mat(row, fromCol + index); }
        [[nodiscard]] __device__ const ScalarType& operator[](size_t index) const { assert(index < getLength()); return mat(row, fromCol + index); }
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Length == Dynamic ? colCount : Length; }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t index) { return mat.data_ptr(row, fromCol + index); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t index) const { return mat.data_ptr(row, fromCol + index); }
    };

    template<class MatrixType, size_t Length>
    class device_obj<ColContinuousVector<MatrixType, Length>> : public Traits<device_obj<ColContinuousVector<MatrixType, Length>>>::Base {
    public:
        using host_obj = ColContinuousVector<MatrixType, Length>;
        using This = device_obj<host_obj>;
        using Base = typename Traits<This>::Base;
        using ScalarType = typename MatrixType::ScalarType;
    private:
        device_obj<MatrixType>& mat;
        size_t col;
        size_t fromRow;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : mat(mat_.getDerived()), col(col_), fromRow(fromRow_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(col < mat.getColumn());
        }
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        __device__ device_obj& operator=(const device_obj& v) { v.assignTo(*this); return *this; }
        __device__ device_obj& operator=(device_obj&& v) noexcept { return operator=(std::cref(v)); }
        [[nodiscard]] __device__ ScalarType& operator[](size_t index) { assert(index < getLength()); return mat(fromRow + index, col); }
        [[nodiscard]] __device__ const ScalarType& operator[](size_t index) const { assert(index < getLength()); return mat(fromRow + index, col); }
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Length == Dynamic ? rowCount : Length; }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t index) { return mat.data_ptr(fromRow + index, col); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t index) const { return mat.data_ptr(fromRow + index, col); }
    };

    template<class MatrixType, size_t Column>
    class device_obj<ContinuousMatrixBlock<MatrixType, 1, Column>>
            : public device_obj<LValueMatrix<ContinuousMatrixBlock<MatrixType, 1, Column>>>
            , public device_obj<RowContinuousVector<MatrixType, Column>> {
        using host_obj = ContinuousMatrixBlock<MatrixType, 1, Column>;
        using This = device_obj<host_obj>;
    public:
        using Base = device_obj<LValueMatrix<host_obj>>;
        using VectorBase = device_obj<RowContinuousVector<MatrixType, Column>>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t fromRow_, [[maybe_unused]] size_t rowCount_, size_t fromCol_, size_t colCount_)
                : device_obj(mat_, fromRow_, fromCol_, colCount_) {
            assert(rowCount_ == 1);
        }
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t row_, size_t fromCol_, size_t colCount_)
                : VectorBase(mat_, row_, fromCol_, colCount_) {}
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        __device__ device_obj& operator=(const device_obj& m) { VectorBase::operator=(m.asVector()); return *this; }
        __device__ device_obj& operator=(device_obj&& m) noexcept { VectorBase::operator=(m.asVector()); return *this; }
        template<class T> This& operator=(const ScalarBase<T>& s) { VectorBase::operator=(s); return *this; }
        using Base::operator=;
        using VectorBase::operator=;
        [[nodiscard]] __device__ ScalarType& operator()([[maybe_unused]] size_t row, size_t col) { assert(row == 0); return VectorBase::operator[](col); }
        [[nodiscard]] __device__ const ScalarType& operator()([[maybe_unused]] size_t row, size_t col) const { assert(row == 0); return VectorBase::operator[](col); }
        template<class T> void operator+=(const ScalarBase<T>& s) { VectorBase::operator+=(s); }
        template<class T> void operator-=(const ScalarBase<T>& s) { VectorBase::operator-=(s); }
        template<class T> void operator*=(const ScalarBase<T>& s) { VectorBase::operator*=(s); }
        template<class T> void operator/=(const ScalarBase<T>& s) { VectorBase::operator/=(s); }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == getColumn()); }

        using Base::row;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        using VectorBase::data_ptr;
        [[nodiscard]] __host__ __device__ constexpr static size_t getRow() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return Column == Dynamic ? VectorBase::getLength() : Column; }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr([[maybe_unused]] size_t row, size_t column) { assert(row == 0); return VectorBase::data_ptr(column); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr([[maybe_unused]] size_t row, size_t column) const { assert(row == 0); return VectorBase::data_ptr(column); }
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
            , public device_obj<ColContinuousVector<MatrixType, Row>> {
        using host_obj = ContinuousMatrixBlock<MatrixType, Row, 1>;
        using This = device_obj<host_obj>;
    public:
        using Base = device_obj<LValueMatrix<host_obj>>;
        using VectorBase = device_obj<ColContinuousVector<MatrixType, Row>>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, [[maybe_unused]] size_t colCount_)
                : device_obj(mat_, fromRow_, rowCount_, fromCol_) {
            assert(colCount_ == 1);
        }
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        __device__ device_obj& operator=(const device_obj& m) { VectorBase::operator=(m.asVector()); return *this; }
        __device__ device_obj& operator=(device_obj&& m) noexcept { VectorBase::operator=(m.asVector()); return *this; }
        template<class T> This& operator=(const ScalarBase<T>& s) { VectorBase::operator=(s); return *this; }
        using Base::operator=;
        using VectorBase::operator=;
        [[nodiscard]] __device__ ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) { assert(col == 0); return VectorBase::operator[](row); }
        [[nodiscard]] __device__ const ScalarType& operator()(size_t row, [[maybe_unused]] size_t col) const { assert(col == 0); return VectorBase::operator[](row); }
        template<class T> void operator+=(const ScalarBase<T>& s) { VectorBase::operator+=(s); }
        template<class T> void operator-=(const ScalarBase<T>& s) { VectorBase::operator-=(s); }
        template<class T> void operator*=(const ScalarBase<T>& s) { VectorBase::operator*=(s); }
        template<class T> void operator/=(const ScalarBase<T>& s) { VectorBase::operator/=(s); }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == getRow() && col == 1); }

        using Base::col;
        /* Getters */
        using Base::calc;
        using VectorBase::calc;
        using VectorBase::max;
        using VectorBase::min;
        using VectorBase::sum;
        using VectorBase::data_ptr;
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Row == Dynamic ? VectorBase::getLength() : Row; }
        [[nodiscard]] __host__ __device__ constexpr static size_t getColumn() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ ScalarType* data_ptr(size_t row, [[maybe_unused]] size_t column) { assert(column == 0); return VectorBase::data_ptr(row); }
        [[nodiscard]] __host__ __device__ const ScalarType* data_ptr(size_t row, [[maybe_unused]] size_t column) const { assert(column == 0); return VectorBase::data_ptr(row); }
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
            , public device_obj<ColContinuousVector<MatrixType, 1>> {
        using host_obj = ContinuousMatrixBlock<MatrixType, 1, 1>;
        using This = device_obj<host_obj>;
    public:
        using Base = device_obj<LValueMatrix<ContinuousMatrixBlock<MatrixType, 1, 1>>>;
        using VectorBase = device_obj<ColContinuousVector<MatrixType, 1>>;
        using ScalarType = typename MatrixType::ScalarType;
    public:
        __host__ __device__ device_obj(device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t col_)
                : VectorBase(mat_, fromRow_, rowCount_, col_) {}
        device_obj(const device_obj&) = delete;
        device_obj(device_obj&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        __device__ device_obj& operator=(const device_obj& m) { VectorBase::operator=(m.asVector()); return *this; }
        __device__ device_obj& operator=(device_obj&& m) noexcept { VectorBase::operator=(m.asVector()); return *this; }
        using Base::operator=;
        using VectorBase::operator=;
        [[nodiscard]] __device__ ScalarType& operator()([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
            assert(row == 0 && col == 0);
            return VectorBase::operator[](0);
        }
        [[nodiscard]] __device__ const ScalarType& operator()([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) const {
            assert(row == 0 && col == 0);
            return VectorBase::operator[](0);
        }
        /* Operations */
        using Base::assignTo;
        using VectorBase::assignTo;
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == 1 && col == 1); }
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
        [[nodiscard]] __host__ __device__ Base& asMatrix() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const Base& asMatrix() const noexcept { return *this; }
        [[nodiscard]] __host__ __device__ VectorBase& asVector() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const VectorBase& asVector() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row, size_t Column>
    class device_obj<ContinuousMatrixBlock<MatrixType, Row, Column>>
            : public device_obj<LValueMatrix<ContinuousMatrixBlock<MatrixType, Row, Column>>> {
        using host_obj = ContinuousMatrixBlock<MatrixType, Row, Column>;
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
        __device__ device_obj& operator=(const device_obj& m) { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        __device__ device_obj& operator=(device_obj&& m) noexcept { Base::operator=(static_cast<const typename Base::Base&>(m)); return *this; }
        [[nodiscard]] __device__ ScalarType& operator()(size_t row, size_t col);
        [[nodiscard]] __device__ const ScalarType& operator()(size_t row, size_t col) const;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return Column == Dynamic ? colCount : Column; }
        [[nodiscard]] __host__ __device__ This& asMatrix() noexcept { return *this; }
        [[nodiscard]] __host__ __device__ const This& asMatrix() const noexcept { return *this; }
    };

    template<class MatrixType, size_t Row, size_t Column>
    __host__ __device__ device_obj<ContinuousMatrixBlock<MatrixType, Row, Column>>::device_obj(
            device_obj<ContinuousMatrix<MatrixType>>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(mat_.getDerived())
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getColumn());
    }

    template<class MatrixType, size_t Row, size_t Column>
    __device__ typename device_obj<ContinuousMatrixBlock<MatrixType, Row, Column>>::ScalarType&
    device_obj<ContinuousMatrixBlock<MatrixType, Row, Column>>::operator()(size_t row, size_t col) {
        assert(row < rowCount);
        assert(col < colCount);
        return mat(row + fromRow, col + fromCol);
    }

    template<class MatrixType, size_t Row, size_t Column>
    __device__ const typename device_obj<ContinuousMatrixBlock<MatrixType, Row, Column>>::ScalarType&
    device_obj<ContinuousMatrixBlock<MatrixType, Row, Column>>::operator()(size_t row, size_t col) const {
        assert(row < rowCount);
        assert(col < colCount);
        return mat(row + fromRow, col + fromCol);
    }
}

namespace Physica {
    using namespace Core;

    template<class MatrixType, size_t Length>
    class Traits<Core::device_obj<RowContinuousVector<MatrixType, Length>>> {
        using VectorType = RowContinuousVector<MatrixType, Length>;
        constexpr static bool isRowMatrix = MatrixOption::isRowMatrix<MatrixType>();
        constexpr static bool isElementMatrix = MatrixOption::isElementMatrix<MatrixType>();
        constexpr static bool isRowVector = isElementMatrix && MatrixType::RowAtCompile == 1;
    public:
        using Base = Core::device_obj<typename std::conditional<isRowMatrix || isRowVector, ContinuousVector<VectorType>, LValueVector<VectorType>>::type>;
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static size_t MaxSizeAtCompile = Length;
    };

    template<class MatrixType, size_t Length>
    class Traits<Core::device_obj<ColContinuousVector<MatrixType, Length>>> {
        using VectorType = ColContinuousVector<MatrixType, Length>;
        constexpr static bool isColumnMatrix = MatrixOption::isColumnMatrix<MatrixType>();
        constexpr static bool isElementMatrix = MatrixOption::isElementMatrix<MatrixType>();
        constexpr static bool isColumnVector = isElementMatrix && MatrixType::ColumnAtCompile == 1;
    public:
        using Base = Core::device_obj<typename std::conditional<isColumnMatrix || isColumnVector, ContinuousVector<VectorType>, LValueVector<VectorType>>::type>;
        using ScalarType = typename MatrixType::ScalarType;
        constexpr static size_t SizeAtCompile = Length;
        constexpr static size_t MaxSizeAtCompile = Length;
    };

    template<class MatrixType, size_t Row, size_t Column>
    class Traits<Core::device_obj<ContinuousMatrixBlock<MatrixType, Row, Column>>> : public Traits<ContinuousMatrixBlock<MatrixType, Row, Column>> {};
}
