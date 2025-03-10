/*
 * Copyright 2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/LValueVector.cuh"
#include "../LValueMatrix.cuh"

namespace Physica {
    template<Matrix T, size_t Col>
    class device_obj<LMatrixBlock<T, 1, Col>> : public device_obj<LValueVector<LMatrixBlock<T, 1, Col>>> {
        using host_obj = LMatrixBlock<T, 1, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        PlainStruct<device_obj<T>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<T>& mat_, size_t fromRow_, size_t fromCol_, size_t colCount_)
                : mat(asStruct(mat_)), fromRow(fromRow_), fromCol(fromCol_), colCount(colCount_) {
            assert(fromRow < mat_.getRow());
            assert(fromCol + colCount <= mat_.getCol());
        }
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < colCount);
            return mat.getDerived().data_ptr(fromRow, fromCol + index);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    };

    template<Matrix T, size_t Row>
    class device_obj<LMatrixBlock<T, Row, 1>> : public device_obj<LValueVector<LMatrixBlock<T, Row, 1>>> {
        using host_obj = LMatrixBlock<T, Row, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        PlainStruct<device_obj<T>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<T>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(asStruct(mat_)), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t index) {
            assert(index < rowCount);
            return mat.getDerived().data_ptr(fromRow + index, fromCol);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    };

    template<Matrix T>
    class device_obj<LMatrixBlock<T, 1, 1>> : public device_obj<LValueVector<LMatrixBlock<T, 1, 1>>> {
        using host_obj = LMatrixBlock<T, 1, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
    public:
        using typename Base::ScalarType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        PlainStruct<device_obj<T>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<T>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
                : mat(asStruct(mat_)), fromRow(fromRow_), fromCol(fromCol_), rowCount(rowCount_) {
            assert(fromRow + rowCount <= mat.getRow());
            assert(fromCol < mat.getCol());
        }
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ PtrTy data_ptr([[maybe_unused]] size_t index) {
            assert(index == 0 && "[Error]: Index overflow");
            return mat.getDerived().data_ptr(fromRow, fromCol);
        }
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr([[maybe_unused]] size_t index) const {
            return const_cast<This&>(*this).data_ptr(index);
        }
    };

    template<Matrix T>
    class device_obj<LMatrixBlock<T, Dynamic, Dynamic>> : public device_obj<LValueMatrix<LMatrixBlock<T, Dynamic, Dynamic>>> {
        using host_obj = LMatrixBlock<T, Dynamic, Dynamic>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
    public:
        using typename Base::ScalarType;
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using PtrTy = ScalarType::PtrTy;
        using ConstPtrTy = ScalarType::ConstPtrTy;
    private:
        PlainStruct<device_obj<T>> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<T>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return colCount; }
        [[nodiscard]] __host__ __device__ inline PtrTy data_ptr(size_t row, size_t col);
        [[nodiscard]] __host__ __device__ inline ConstPtrTy data_ptr(size_t row, size_t col) const;
    };

    template<Matrix T>
    __host__ __device__ device_obj<LMatrixBlock<T, Dynamic, Dynamic>>::device_obj(
            device_obj<T>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(asStruct(mat_))
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<Matrix T>
    __host__ __device__ inline auto device_obj<LMatrixBlock<T, Dynamic, Dynamic>>::data_ptr(size_t row, size_t col) -> PtrTy {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.getDerived().data_ptr(row + fromRow, col + fromCol);
    }

    template<Matrix T>
    __host__ __device__ auto device_obj<LMatrixBlock<T, Dynamic, Dynamic>>::data_ptr(size_t row, size_t col) const -> ConstPtrTy {
        return const_cast<This&>(*this).data_ptr(row, col);
    }
}

namespace Physica {
    template<Matrix T, size_t Row, size_t Col>
    class Traits<device_obj<LMatrixBlock<T, Row, Col>>> : public Traits<LMatrixBlock<T, Row, Col>> {};
}
