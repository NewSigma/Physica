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

#include "ContinuousMatrixBlock.h"

namespace Physica {
    template<Matrix M, size_t Col>
    class device_obj<ContinuousMatrixBlock<M, 1, Col>> : public device_obj<ContinuousVector<ContinuousMatrixBlock<M, 1, Col>>> {
        using host_obj = ContinuousMatrixBlock<M, 1, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<ContinuousVector<host_obj>>;
    private:
        PlainStruct<device_obj<M>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat, size_t fromRow, [[maybe_unused]] size_t rowCount, size_t fromCol, size_t colCount);
        __host__ __device__ device_obj(device_obj<M>& mat, size_t fromRow, size_t fromCol, size_t colCount);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ auto data(this auto&& self) noexcept;
    };

    template<Matrix M, size_t Col>
    __host__ __device__ device_obj<ContinuousMatrixBlock<M, 1, Col>>::device_obj(device_obj<M>& mat, size_t fromRow, [[maybe_unused]] size_t rowCount, size_t fromCol, size_t colCount)
            : device_obj(mat, fromRow, fromCol, colCount) {
        assert(rowCount == 1);
    }

    template<Matrix M, size_t Col>
    __host__ __device__ device_obj<ContinuousMatrixBlock<M, 1, Col>>::device_obj(device_obj<M>& mat, size_t fromRow, size_t fromCol, size_t colCount)
            : mat(asStruct(mat)), fromRow(fromRow), fromCol(fromCol), colCount(colCount) {
        assert(fromRow < mat.getRow());
        assert(fromCol + colCount <= mat.getCol());
    }

    template<Matrix M, size_t Col>
    __host__ __device__ auto device_obj<ContinuousMatrixBlock<M, 1, Col>>::data(this auto&& self) noexcept {
        return self.mat.getDerived().data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row>
    class device_obj<ContinuousMatrixBlock<M, Row, 1>> : public device_obj<ContinuousVector<ContinuousMatrixBlock<M, Row, 1>>> {
        using host_obj = ContinuousMatrixBlock<M, Row, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<ContinuousVector<host_obj>>;
    private:
        PlainStruct<device_obj<M>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat, size_t fromRow, size_t rowCount, size_t fromCol, [[maybe_unused]] size_t colCount);
        __host__ __device__ device_obj(device_obj<M>& mat, size_t fromRow, size_t rowCount, size_t fromCol);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getLength() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ auto data(this auto&& self) noexcept;
    };

    template<Matrix M, size_t Row>
    __host__ __device__ device_obj<ContinuousMatrixBlock<M, Row, 1>>::device_obj(device_obj<M>& mat, size_t fromRow, size_t rowCount, size_t fromCol, [[maybe_unused]] size_t colCount)
            : device_obj(mat, fromRow, rowCount, fromCol) {
        assert(colCount == 1);
    }

    template<Matrix M, size_t Row>
    __host__ __device__ device_obj<ContinuousMatrixBlock<M, Row, 1>>::device_obj(device_obj<M>& mat, size_t fromRow, size_t rowCount, size_t fromCol)
            : mat(asStruct(mat)), fromRow(fromRow), fromCol(fromCol), rowCount(rowCount) {
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M, size_t Row>
    __host__ __device__ auto device_obj<ContinuousMatrixBlock<M, Row, 1>>::data(this auto&& self) noexcept {
        return self.mat.getDerived().data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M>
    class device_obj<ContinuousMatrixBlock<M, 1, 1>> : public device_obj<ContinuousVector<ContinuousMatrixBlock<M, 1, 1>>> {
        using host_obj = ContinuousMatrixBlock<M, 1, 1>;
        using This = device_obj<host_obj>;
        using Base = device_obj<ContinuousVector<host_obj>>;
    private:
        PlainStruct<device_obj<M>> mat;
        size_t fromRow;
        size_t fromCol;
        size_t rowCount;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat, size_t fromRow, size_t rowCount, size_t fromCol);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t length) { assert(length == getLength()); }
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr static size_t getLength() noexcept { return 1; }
        [[nodiscard]] __host__ __device__ auto data(this auto&& self) noexcept;
    };

    template<Matrix M>
    __host__ __device__ device_obj<ContinuousMatrixBlock<M, 1, 1>>::device_obj(device_obj<M>& mat, size_t fromRow, size_t rowCount, size_t fromCol)
            : mat(asStruct(mat)), fromRow(fromRow), fromCol(fromCol), rowCount(rowCount) {
        assert(fromRow + rowCount <= mat.getRow());
        assert(fromCol < mat.getCol());
    }

    template<Matrix M>
    __host__ __device__ auto device_obj<ContinuousMatrixBlock<M, 1, 1>>::data(this auto&& self) noexcept {
        return self.mat.getDerived().data_ptr(self.fromRow, self.fromCol);
    }

    template<Matrix M, size_t Row, size_t Col>
    class device_obj<ContinuousMatrixBlock<M, Row, Col>>
            : public device_obj<LValueMatrix<ContinuousMatrixBlock<M, Row, Col>>> {
        using host_obj = ContinuousMatrixBlock<M, Row, Col>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<host_obj>>;
    private:
        PlainStruct<device_obj<M>> mat;
        size_t fromRow;
        size_t rowCount;
        size_t fromCol;
        size_t colCount;
    public:
        __host__ __device__ device_obj(device_obj<M>& mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        device_obj(const This&) = delete;
        device_obj(This&&) noexcept = delete;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& m) { Base::operator=(m); return *this; }
        This& operator=(This&& m) { return *this = m; }
        using Base::operator=;
        /* Operations */
        using Base::resize;
        __host__ __device__ void resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) { assert(row == rowCount && col == colCount); }
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Row == Dynamic ? rowCount : Row; }
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Col == Dynamic ? colCount : Col; }
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&& self, size_t row, size_t col) noexcept;
    };

    template<Matrix M, size_t Row, size_t Col>
    __host__ __device__ device_obj<ContinuousMatrixBlock<M, Row, Col>>::device_obj(device_obj<M>& mat, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount)
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
    __host__ __device__ auto device_obj<ContinuousMatrixBlock<M, Row, Col>>::data_ptr(this auto&& self, size_t row, size_t col) noexcept {
        assert(row < self.rowCount);
        assert(col < self.colCount);
        return self.mat.getDerived().data_ptr(row + self.fromRow, col + self.fromCol);
    }
}

namespace Physica {
    template<Matrix M, size_t Row, size_t Col>
    class Traits<device_obj<ContinuousMatrixBlock<M, Row, Col>>> : public Traits<ContinuousMatrixBlock<M, Row, Col>> {};
}
