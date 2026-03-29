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
        using Ref = add_device_obj<M>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(Ref mat_, size_t fromRow_, size_t fromCol_, size_t colCount_)
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
        using Ref = add_device_obj<M>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(Ref mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_)
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
        using Ref = add_device_obj<M>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        size_t fromCol;
    public:
        __host__ __device__ device_obj(Ref mat_, size_t fromRow_, size_t fromCol_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == 1); }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, [[maybe_unused]] size_t index) noexcept;
    };

    template<Matrix M>
    __host__ __device__ device_obj<LMatrixBlock<M, 1, 1>>::device_obj(Ref mat_, size_t fromRow_, size_t fromCol_)
            : mat(asStruct(mat_)), fromRow(fromRow_), fromCol(fromCol_) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
    }

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
        using Ref = add_device_obj<M>::type;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(Ref mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }

        [[nodiscard]] __host__ __device__ auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] __host__ __device__ auto col(this auto&&, size_t c) noexcept;
        [[nodiscard]] __host__ __device__ auto rows(this auto&&, size_t fromRow, size_t rowCount) noexcept;
        [[nodiscard]] __host__ __device__ auto topRows(this auto&&, size_t to) noexcept;
        [[nodiscard]] __host__ __device__ auto bottomRows(this auto&&, size_t from) noexcept;
        [[nodiscard]] __host__ __device__ auto cols(this auto&&, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] __host__ __device__ auto leftCols(this auto&&, size_t to) noexcept;
        [[nodiscard]] __host__ __device__ auto rightCols(this auto&&, size_t from) noexcept;
        [[nodiscard]] __host__ __device__ auto topLeftCorner(this auto&&, size_t toRow, size_t toCol) noexcept;
        [[nodiscard]] __host__ __device__ auto topLeftCorner(this auto&&, size_t to) noexcept;
        [[nodiscard]] __host__ __device__ auto topRightCorner(this auto&&, size_t toRow, size_t fromCol) noexcept;
        [[nodiscard]] __host__ __device__ auto bottomLeftCorner(this auto&&, size_t fromRow, size_t toCol) noexcept;
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(this auto&&, size_t fromRow, size_t fromCol) noexcept;
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(this auto&&, size_t from) noexcept;
        [[nodiscard]] __host__ __device__ auto block(this auto&&, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return rowCount; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return colCount; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
    };

    template<Matrix M>
    __host__ __device__ device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::device_obj(
            Ref mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(asStruct(mat_))
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::row(this auto&& self, size_t r) noexcept {
        assert(r < self.getRow());
        return device_obj<LMatrixBlock<M, 1, Dynamic>>(self.mat, self.fromRow + r, self.fromCol, self.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::col(this auto&& self, size_t c) noexcept {
        assert(c < self.getCol());
        return device_obj<LMatrixBlock<M, Dynamic, 1>>(self.mat, self.fromRow, self.getRow(), self.fromCol + c);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        check(self, fromRow, rowCount, 0, self.getCol());
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + fromRow, rowCount, self.fromCol, self.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::topRows(this auto&& self, size_t to) noexcept {
        check(self, 0, to, 0, self.getCol());
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, to, self.fromCol, self.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::bottomRows(this auto&& self, size_t from) noexcept {
        check(self, from, self.getRow() - from, 0, self.getCol());
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + from, self.getRow() - from, self.fromCol, self.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        check(self, 0, self.getRow(), fromCol, colCount);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, self.getRow(), self.fromCol + fromCol, colCount);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::leftCols(this auto&& self, size_t to) noexcept {
        check(self, 0, self.getRow(), 0, to);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, self.getRow(), self.fromCol, to);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::rightCols(this auto&& self, size_t from) noexcept {
        check(self, 0, self.getRow(), from, self.getCol() - from);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, self.getRow(), self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        check(self, 0, toRow, 0, toCol);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, toRow, self.fromCol, toCol);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::topLeftCorner(this auto&& self, size_t to) noexcept {
        check(self, 0, to, 0, to);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, to, self.fromCol, to);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        check(self, 0, toRow, fromCol, self.getCol() - fromCol);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, toRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        check(self, fromRow, self.getRow() - fromRow, 0, toCol);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol, toCol);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        check(self, fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        check(self, from, self.getRow() - from, from, self.getCol() - from);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + from, self.getRow() - from, self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<LMatrixBlock<M, Dynamic, Dynamic>>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        check(self, fromRow, rowCount, fromCol, colCount);
        return device_obj<LMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + fromRow, rowCount, self.fromCol + fromCol, colCount);
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
    class Traits<device_obj<LMatrixBlock<M, Row, Col>>> : public Traits<LMatrixBlock<M, Row, Col>> {
        static_assert(!DeviceObj<M>);
    };
}
