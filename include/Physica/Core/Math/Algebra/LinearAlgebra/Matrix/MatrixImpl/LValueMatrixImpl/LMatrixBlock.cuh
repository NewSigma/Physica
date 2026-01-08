/*
 * Copyright 2025-2026 Weibo He.
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
    template<Matrix M, size_t Col>
    class device_obj<LMatrixBlock<M, 1, Col>> : public device_obj<LValueVector<LMatrixBlock<M, 1, Col>>> {
        using host_obj = LMatrixBlock<M, 1, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
    private:
        PlainStruct<device_obj<M>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat_, size_t fromRow_, size_t fromCol_, size_t colCount_)
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
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t index) noexcept;
    };

    template<Matrix M, size_t Col>
    __host__ __device__ auto device_obj<LMatrixBlock<M, 1, Col>>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.colCount);
        return self.mat.getDerived().data_ptr(self.fromRow, self.fromCol + index);
    }

    template<Matrix M, size_t Row>
    class device_obj<LMatrixBlock<M, Row, 1>> : public device_obj<LValueVector<LMatrixBlock<M, Row, 1>>> {
        using host_obj = LMatrixBlock<M, Row, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
    private:
        PlainStruct<device_obj<M>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
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
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t index) noexcept;
    };

    template<Matrix M, size_t Row>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Row, 1>>::data_ptr(this auto&& self, size_t index) noexcept {
        assert(index < self.rowCount);
        return self.mat.getDerived().data_ptr(self.fromRow + index, self.fromCol);
    }

    template<Matrix M>
    class device_obj<LMatrixBlock<M, 1, 1>> : public device_obj<LValueVector<LMatrixBlock<M, 1, 1>>> {
        using host_obj = LMatrixBlock<M, 1, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueVector<host_obj>>;
    private:
        PlainStruct<device_obj<M>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
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
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == rowCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, [[maybe_unused]] size_t index) noexcept;
    };

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, 1, 1>>::data_ptr(this auto&& self, [[maybe_unused]] size_t index) noexcept {
        assert(index == 0 && "[Error]: Index overflow");
        return self.mat.getDerived().data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M>
    class device_obj<LMatrixBlock<M, Dynamic, Dynamic>> : public device_obj<LValueMatrix<LMatrixBlock<M, Dynamic, Dynamic>>> {
        using host_obj = LMatrixBlock<M, Dynamic, Dynamic>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
    private:
        PlainStruct<device_obj<M>> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return colCount; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
    };

    template<Matrix M>
    __host__ __device__ device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::device_obj(
            device_obj<M>& mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(asStruct(mat_))
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.rowCount);
        assert(col < self.colCount);
        return self.mat.getDerived().data_ptr(row + self.fromRow, col + self.fromCol);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<device_obj<LMatrixBlock<M, Row, Col>>> : public Traits<LMatrixBlock<M, Row, Col>> {};
}
