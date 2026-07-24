/*
 * Copyright 2023-2026 Weibo He.
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

#include "CompactMatrixBlock.h"
#include "../CompactMatrix.cuh"
#include "Physica/Core/Utils/Empty.h"

namespace Physica {
    template<Matrix M, size_t Row, size_t Col> requires (Row == 1 || Col == 1)
    class device_obj<CompactMatrixBlock<M, Row, Col>> : public device_obj<CompactVector<CompactMatrixBlock<M, Row, Col>>> {
        static_assert((Row == 1 && MatrixMajor::isRowMatrix<M>()) || (Col == 1 && MatrixMajor::isColMatrix<M>()) || (Row == 1 || Col == 1),
                      "[Error]: The matrix block is not compact");
        using host_obj = CompactMatrixBlock<M, Row, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<CompactVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
        using MaybeRowCount = std::conditional_t<Row == Dynamic, size_t, Empty>;
        using MaybeColCount = std::conditional_t<Col == Dynamic, size_t, Empty>;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        [[no_unique_address]] MaybeRowCount rowCount;
        size_t fromCol;
        [[no_unique_address]] MaybeColCount colCount;
    public:
        __host__ __device__ device_obj(Ref mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& m);
        This& operator=(This&& m);
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize(size_t length);
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept;
        [[nodiscard]] __host__ __device__ auto data(this auto&& self) noexcept;
    };

    template<Matrix M, size_t Row, size_t Col> requires (Row == 1 || Col == 1)
    __host__ __device__ device_obj<CompactMatrixBlock<M, Row, Col>>::device_obj(Ref mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
            : mat(asStruct(mat))
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol)
            , colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol + colCount <= mat.getCol());
        if constexpr (Row != Dynamic)
            assert(Row == rowCount);
        if constexpr (Col != Dynamic)
            assert(Col == colCount);
    }

    template<Matrix M, size_t Row, size_t Col> requires (Row == 1 || Col == 1)
    auto device_obj<CompactMatrixBlock<M, Row, Col>>::operator=(const This& m) -> This& {
        Base::operator=(m);
        return *this;
    }

    template<Matrix M, size_t Row, size_t Col> requires (Row == 1 || Col == 1)
    auto device_obj<CompactMatrixBlock<M, Row, Col>>::operator=(This&& m) -> This& {
        return *this = m;
    }

    template<Matrix M, size_t Row, size_t Col> requires (Row == 1 || Col == 1)
    __host__ __device__ void device_obj<CompactMatrixBlock<M, Row, Col>>::resize([[maybe_unused]] size_t length) {
        assert(length == getLength());
    }

    template<Matrix M, size_t Row, size_t Col> requires (Row == 1 || Col == 1)
    __host__ __device__ size_t device_obj<CompactMatrixBlock<M, Row, Col>>::getLength() const noexcept {
        if constexpr (Row == 1) {
            if constexpr (Col == Dynamic)
                return colCount;
            else
                return Col;
        }
        else {
            if constexpr (Row == Dynamic)
                return rowCount;
            else
                return Row;
        }
    }

    template<Matrix M, size_t Row, size_t Col> requires (Row == 1 || Col == 1)
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::data(this auto&& self) noexcept {
        return self.mat.getDerived().data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    class device_obj<CompactMatrixBlock<M, Row, Col>>
            : public device_obj<LValueMatrix<CompactMatrixBlock<M, Row, Col>>> {
        using host_obj = CompactMatrixBlock<M, Row, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
        using Ref = add_device_obj<M>::type;
        using MaybeRowCount = std::conditional_t<Row == Dynamic, size_t, Empty>;
        using MaybeColCount = std::conditional_t<Col == Dynamic, size_t, Empty>;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        [[no_unique_address]] MaybeRowCount rowCount;
        size_t fromCol;
        [[no_unique_address]] MaybeColCount colCount;
    public:
        __host__ __device__ device_obj(Ref mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& m);
        This& operator=(This&& m);
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize(size_t row, size_t col);

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
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept;
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t row, size_t col) noexcept;
    };

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ device_obj<CompactMatrixBlock<M, Row, Col>>::device_obj(Ref mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
            : mat(asStruct(mat))
            , fromRow(fromRow)
            , rowCount(rowCount)
            , fromCol(fromCol)
            , colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol < mat.getCol());
        assert((fromRow + rowCount) <= mat.getRow());
        assert((fromCol + colCount) <= mat.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    auto device_obj<CompactMatrixBlock<M, Row, Col>>::operator=(const This& m) -> This& {
        Base::operator=(m);
        return *this;
    }

    template<Matrix M, size_t Row, size_t Col>
    auto device_obj<CompactMatrixBlock<M, Row, Col>>::operator=(This&& m) -> This& {
        return *this = m;
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ void device_obj<CompactMatrixBlock<M, Row, Col>>::resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::row(this auto&& self, size_t r) noexcept {
        assert(r < self.getRow());
        if constexpr (MatrixMajor::isRowMatrix<M>())
            return device_obj<CompactMatrixBlock<M, 1, Col>>(self.mat, self.fromRow + r, 1, self.fromCol, self.getCol());
        else
            return device_obj<LMatrixBlock<M, 1, Col>>(self.mat, self.fromRow + r, 1, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::col(this auto&& self, size_t c) noexcept {
        assert(c < self.getCol());
        if constexpr (MatrixMajor::isColMatrix<M>())
            return device_obj<CompactMatrixBlock<M, Row, 1>>(self.mat, self.fromRow, self.getRow(), self.fromCol + c, 1);
        else
            return device_obj<LMatrixBlock<M, Row, 1>>(self.mat, self.fromRow, self.getRow(), self.fromCol + c, 1);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        check(self, fromRow, rowCount, 0, self.getCol());
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow + fromRow, rowCount, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::topRows(this auto&& self, size_t to) noexcept {
        check(self, 0, to, 0, self.getCol());
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow, to, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::bottomRows(this auto&& self, size_t from) noexcept {
        check(self, from, self.getRow() - from, 0, self.getCol());
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow + from, self.getRow() - from, self.fromCol, self.getCol());
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        check(self, 0, self.getRow(), fromCol, colCount);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow, self.getRow(), self.fromCol + fromCol, colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::leftCols(this auto&& self, size_t to) noexcept {
        check(self, 0, self.getRow(), 0, to);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow, self.getRow(), self.fromCol, to);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::rightCols(this auto&& self, size_t from) noexcept {
        check(self, 0, self.getRow(), from, self.getCol() - from);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow, self.getRow(), self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        check(self, 0, toRow, 0, toCol);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow, toRow, self.fromCol, toCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::topLeftCorner(this auto&& self, size_t to) noexcept {
        check(self, 0, to, 0, to);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow, to, self.fromCol, to);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        check(self, 0, toRow, fromCol, self.getCol() - fromCol);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow, toRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        check(self, fromRow, self.getRow() - fromRow, 0, toCol);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol, toCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        check(self, fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        check(self, from, self.getRow() - from, from, self.getCol() - from);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow + from, self.getRow() - from, self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        check(self, fromRow, rowCount, fromCol, colCount);
        return device_obj<CompactMatrixBlock<M, Row, Col>>(self.mat, self.fromRow + fromRow, rowCount, self.fromCol + fromCol, colCount);
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ size_t device_obj<CompactMatrixBlock<M, Row, Col>>::getRow() const noexcept {
        if constexpr (Row == Dynamic)
            return rowCount;
        else
            return Row;
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ size_t device_obj<CompactMatrixBlock<M, Row, Col>>::getCol() const noexcept {
        if constexpr (Col == Dynamic)
            return colCount;
        else
            return Col;
    }

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ auto device_obj<CompactMatrixBlock<M, Row, Col>>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.getRow());
        assert(col < self.getCol());
        return self.mat.getDerived().data_ptr(row + self.fromRow, col + self.fromCol);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<device_obj<CompactMatrixBlock<M, Row, Col>>> : public Traits<CompactMatrixBlock<M, Row, Col>> {
        static_assert(!is_device_obj<M>::value);
    };
}
