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

namespace Physica {
    template<class Derived>
    class device_obj<ContinuousMatrix<Derived>> : public device_obj<LValueMatrix<Derived>> {
        using host_obj = ContinuousMatrix<Derived>;
        using This = device_obj<host_obj>;
        using Base = device_obj<LValueMatrix<Derived>>;
    public:
        using typename Base::ScalarType;
        using Base::RowAtCompile;
        using Base::ColAtCompile;
    protected:
        using typename Base::PtrTy;
        using typename Base::ConstPtrTy;
    private:
        using RowVector = device_obj<ContinuousMatrixBlock<Derived, 1, ColAtCompile>>;
        using ColVector = device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, 1>>;
        template<size_t Row>
        using RowBlock = device_obj<ContinuousMatrixBlock<Derived, Row, ColAtCompile>>;
        template<size_t Col>
        using ColBlock = device_obj<ContinuousMatrixBlock<Derived, RowAtCompile, Col>>;
        template<size_t Row, size_t Col>
        using BlockType = device_obj<ContinuousMatrixBlock<Derived, Row, Col>>;
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Operations */
        [[nodiscard]] __host__ __device__ auto row(size_t r) noexcept;
        [[nodiscard]] __host__ __device__ const auto row(size_t r) const noexcept;
        [[nodiscard]] __host__ __device__ auto col(size_t c) noexcept;
        [[nodiscard]] __host__ __device__ const auto col(size_t c) const noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ auto rows(size_t fromRow, size_t rowCount) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ const auto rows(size_t fromRow, size_t rowCount) const noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ auto topRows(size_t to) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ const auto topRows(size_t to) const noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ auto bottomRows(size_t from) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ const auto bottomRows(size_t from) const noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto cols(size_t fromCol, size_t colCount) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto cols(size_t fromCol, size_t colCount) const noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto leftCols(size_t to) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto leftCols(size_t to) const noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto rightCols(size_t from) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto rightCols(size_t from) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto topLeftCorner(size_t toRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto topLeftCorner(size_t toRow, size_t toCol) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto topLeftCorner(size_t to) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto topLeftCorner(size_t to) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto topRightCorner(size_t toRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto topRightCorner(size_t toRow, size_t fromCol) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto bottomLeftCorner(size_t fromRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(size_t fromRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(size_t from) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto bottomRightCorner(size_t from) const noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ const auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept;

        [[nodiscard]] auto flatten();
        [[nodiscard]] const auto flatten() const;
        /* Getters */
        [[nodiscard]] __host__ __device__ PtrTy data() { return Base::getDerived().data_ptr(0, 0); }
        [[nodiscard]] __host__ __device__ ConstPtrTy data() const { return Base::getDerived().data_ptr(0, 0); }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "ContinuousMatrixImpl/ContinuousMatrixImpl.cuh"
#include "ContinuousMatrixImpl/Flatten.cuh"
