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

#include "RMatrixBlock.h"

namespace Physica {
    template<Matrix M>
    class device_obj<RMatrixBlock<M, 1, Dynamic>> : public device_obj<RValueVector<RMatrixBlock<M, 1, Dynamic>>> {
        using host_obj = RMatrixBlock<M, 1, Dynamic>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
    public:
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(Ref mat_, size_t fromRow, size_t fromCol, size_t colCount);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t index) const {
            assert(index < colCount);
            return mat.getDerived().calc(fromRow, fromCol + index);
        }
        [[nodiscard]] __device__ Tv calc_value(size_t index) const {
            assert(index < colCount);
            return mat.getDerived().calc_value(fromRow, fromCol + index);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return colCount; }
    };

    template<Matrix M>
    __host__ __device__ device_obj<RMatrixBlock<M, 1, Dynamic>>::device_obj(Ref mat_, size_t fromRow, size_t fromCol, size_t colCount)
            : mat(asStruct(mat_)), fromRow(fromRow), fromCol(fromCol), colCount(colCount) {
        assert(fromRow < mat_.getRow());
        assert(fromCol + colCount <= mat_.getCol());
    }

    template<Matrix M>
    class device_obj<RMatrixBlock<M, Dynamic, 1>> : public device_obj<RValueVector<RMatrixBlock<M, Dynamic, 1>>> {
        using host_obj = RMatrixBlock<M, Dynamic, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueVector<host_obj>>;
        using Ref = add_device_obj<M>::type;
    public:
        using Base::isComplex;
        using Base::SizeAtCompile;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(Ref mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Getters */
        [[nodiscard]] __device__ T calc(size_t index) const {
            assert(index < rowCount);
            return mat.getDerived().calc(fromRow + index, fromCol);
        }
        [[nodiscard]] __device__ Tv calc_value(size_t index) const {
            assert(index < rowCount);
            return mat.getDerived().calc_value(fromRow + index, fromCol);
        }
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return rowCount; }
    };

    template<Matrix M>
    __host__ __device__ device_obj<RMatrixBlock<M, Dynamic, 1>>::device_obj(Ref mat_, size_t fromRow, size_t rowCount, size_t fromCol)
            : mat(asStruct(mat_)), fromRow(fromRow), fromCol(fromCol), rowCount(rowCount) {
        assert(fromRow + rowCount <= mat_.getRow());
        assert(fromCol < mat_.getCol());
    }

    template<Matrix M>
    class device_obj<RMatrixBlock<M, Dynamic, Dynamic>> : public device_obj<RValueMatrix<RMatrixBlock<M, Dynamic, Dynamic>>> {
        using host_obj = RMatrixBlock<M, Dynamic, Dynamic>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref = add_device_obj<M>::type;
    protected:
        using typename Base::T;
        using typename Base::Tv;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M>>> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(Ref mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operations */
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const;

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
    };

    template<Matrix M>
    __host__ __device__ device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::device_obj(
            Ref mat_, size_t fromRow_, size_t rowCount_, size_t fromCol_, size_t colCount_)
            : mat(asStruct(mat_))
            , fromRow(fromRow_)
            , rowCount(rowCount_)
            , fromCol(fromCol_)
            , colCount(colCount_) {
        assert((fromRow + rowCount) <= mat_.getRow());
        assert((fromCol + colCount) <= mat_.getCol());
    }

    template<Matrix M>
    __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::calc(size_t row, size_t col) const -> T {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.getDerived().calc(row + fromRow, col + fromCol);
    }

    template<Matrix M>
    __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::calc_value(size_t row, size_t col) const -> Tv {
        assert(row < rowCount);
        assert(col < colCount);
        return mat.getDerived().calc_value(row + fromRow, col + fromCol);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::row(this auto&& self, size_t r) noexcept {
        assert(r < self.getRow());
        return device_obj<RMatrixBlock<M, 1, Dynamic>>(self.mat, self.fromRow + r, self.fromCol, self.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::col(this auto&& self, size_t c) noexcept {
        assert(c < self.getCol());
        return device_obj<RMatrixBlock<M, Dynamic, 1>>(self.mat, self.fromRow, self.getRow(), self.fromCol + c);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::rows(this auto&& self, size_t fromRow, size_t rowCount) noexcept {
        check(self, fromRow, rowCount, 0, self.getCol());
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + fromRow, rowCount, self.fromCol, self.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::topRows(this auto&& self, size_t to) noexcept {
        check(self, 0, to, 0, self.getCol());
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, to, self.fromCol, self.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::bottomRows(this auto&& self, size_t from) noexcept {
        check(self, from, self.getRow() - from, 0, self.getCol());
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + from, self.getRow() - from, self.fromCol, self.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::cols(this auto&& self, size_t fromCol, size_t colCount) noexcept {
        check(self, 0, self.getRow(), fromCol, colCount);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, self.getRow(), self.fromCol + fromCol, colCount);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::leftCols(this auto&& self, size_t to) noexcept {
        check(self, 0, self.getRow(), 0, to);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, self.getRow(), self.fromCol, to);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::rightCols(this auto&& self, size_t from) noexcept {
        check(self, 0, self.getRow(), from, self.getCol() - from);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, self.getRow(), self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::topLeftCorner(this auto&& self, size_t toRow, size_t toCol) noexcept {
        check(self, 0, toRow, 0, toCol);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, toRow, self.fromCol, toCol);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::topLeftCorner(this auto&& self, size_t to) noexcept {
        check(self, 0, to, 0, to);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, to, self.fromCol, to);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::topRightCorner(this auto&& self, size_t toRow, size_t fromCol) noexcept {
        check(self, 0, toRow, fromCol, self.getCol() - fromCol);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow, toRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::bottomLeftCorner(this auto&& self, size_t fromRow, size_t toCol) noexcept {
        check(self, fromRow, self.getRow() - fromRow, 0, toCol);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol, toCol);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::bottomRightCorner(this auto&& self, size_t fromRow, size_t fromCol) noexcept {
        check(self, fromRow, self.getRow() - fromRow, fromCol, self.getCol() - fromCol);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + fromRow, self.getRow() - fromRow, self.fromCol + fromCol, self.getCol() - fromCol);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::bottomRightCorner(this auto&& self, size_t from) noexcept {
        check(self, from, self.getRow() - from, from, self.getCol() - from);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + from, self.getRow() - from, self.fromCol + from, self.getCol() - from);
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<RMatrixBlock<M, Dynamic, Dynamic>>::block(this auto&& self, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept {
        check(self, fromRow, rowCount, fromCol, colCount);
        return device_obj<RMatrixBlock<M, Dynamic, Dynamic>>(self.mat, self.fromRow + fromRow, rowCount, self.fromCol + fromCol, colCount);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<device_obj<RMatrixBlock<M, Row, Col>>> : public Traits<RMatrixBlock<M, Row, Col>> {
        static_assert(!is_device_obj<M>::value);
    };
}
