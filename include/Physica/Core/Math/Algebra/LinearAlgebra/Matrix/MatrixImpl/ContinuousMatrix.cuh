/*
 * Copyright 2023 Weibo He.
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

#include "ContinuousMatrix.h"
#include "LValueMatrix.cuh"
#include "ContinuousMatrixImpl/ContinuousMatrixBlock.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<ContinuousMatrix<Derived>> : public device_obj<LValueMatrix<Derived>> {
        using host_obj = ContinuousMatrix<Derived>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<Derived>>;
    public:
        using typename Base::ScalarType;
        using Base::RowAtCompile;
        using Base::ColumnAtCompile;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    public:
        ~device_obj() = default;
        /* Operators */
        inline device_obj& operator=(const device_obj& obj);
        inline device_obj& operator=(device_obj&& obj) noexcept;
        using Base::operator=;
        /* Operations */
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, 1, ColumnAtCompile>> row(size_t r);
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, 1, ColumnAtCompile>> row(size_t r) const;
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, 1>> col(size_t c);
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, 1>> col(size_t c) const;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, ColumnAtCompile>> rows(size_t fromRow, size_t rowCount);
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, ColumnAtCompile>> rows(size_t fromRow, size_t rowCount) const;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, ColumnAtCompile>> topRows(size_t to);
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, ColumnAtCompile>> topRows(size_t to) const;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, ColumnAtCompile>> bottomRows(size_t from);
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, ColumnAtCompile>> bottomRows(size_t from) const;
        template<size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, Column>> cols(size_t fromCol, size_t colCount);
        template<size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, Column>> cols(size_t fromCol, size_t colCount) const;
        template<size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, Column>> leftCols(size_t to);
        template<size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, Column>> leftCols(size_t to) const;
        template<size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, Column>> rightCols(size_t from);
        template<size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, Column>> rightCols(size_t from) const;
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>> topLeftCorner(size_t toRow, size_t toCol);
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>> topLeftCorner(size_t toRow, size_t toCol) const;
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>> topLeftCorner(size_t to);
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>> topLeftCorner(size_t to) const;
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>> topRightCorner(size_t toRow, size_t fromCol);
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>> topRightCorner(size_t toRow, size_t fromCol) const;
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>> bottomLeftCorner(size_t fromRow, size_t toCol);
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>> bottomLeftCorner(size_t fromRow, size_t toCol) const;
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>> bottomRightCorner(size_t fromRow, size_t fromCol);
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>> bottomRightCorner(size_t fromRow, size_t fromCol) const;
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>> bottomRightCorner(size_t from);
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>> bottomRightCorner(size_t from) const;
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline device_obj<ContinuousMatrixBlock<Derived, Row, Column>> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        template<size_t Row = Dynamic, size_t Column = Dynamic>
        [[nodiscard]] __host__ __device__ inline const device_obj<ContinuousMatrixBlock<Derived, Row, Column>> block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ PtrTy data() { return Base::getDerived().data_ptr(0, 0); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data() const { return Base::getDerived().data_ptr(0, 0); }
        [[nodiscard]] device_obj<ContinuousFlatten<Derived>> flatten() { return {*this}; }
        [[nodiscard]] const device_obj<ContinuousFlatten<Derived>> flatten() const { return {const_cast<This&>(*this)}; }
    protected:
        device_obj() = default;
        device_obj(const device_obj&) = default;
        device_obj(device_obj&&) noexcept = default;
    };
}

#include "ContinuousMatrixImpl/ContinuousMatrixImpl.cuh"
#include "ContinuousMatrixImpl/ContinuousFlatten.cuh"
