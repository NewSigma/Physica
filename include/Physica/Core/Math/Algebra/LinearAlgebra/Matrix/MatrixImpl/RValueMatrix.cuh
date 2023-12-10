/*
 * Copyright 2022-2023 WeiBo He.
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

#include "RMatrixBlock.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<RValueMatrix<Derived>> : public Utils::CRTPBase<device_obj<Derived>> {
        static_assert(!Utils::is_device_obj<Derived>::value, "[Error]: Nested device_obj is unnecessary");
        using Base = Utils::CRTPBase<device_obj<Derived>>;
        using host_obj = RValueMatrix<Derived>;
        using BlockType = device_obj<RMatrixBlock<Derived>>;
    public:
        using ScalarType = typename host_obj::ScalarType;
        using RowVector = device_obj<RMatrixBlock<Derived, 1, Dynamic>>;
        using ColVector = device_obj<RMatrixBlock<Derived, Dynamic, 1>>;
        constexpr static int Option = host_obj::Option;
        constexpr static size_t RowAtCompile = host_obj::RowAtCompile;
        constexpr static size_t ColumnAtCompile = host_obj::ColumnAtCompile;
        constexpr static size_t MaxRowAtCompile = host_obj::MaxRowAtCompile;
        constexpr static size_t MaxColumnAtCompile = host_obj::MaxColumnAtCompile;
        constexpr static size_t SizeAtCompile = host_obj::SizeAtCompile;
        constexpr static size_t MaxSizeAtCompile = host_obj::MaxSizeAtCompile;

        constexpr static bool isColumnMatrix = host_obj::isColumnMatrix;
        constexpr static bool isRowMatrix = host_obj::isRowMatrix;
    public:
        /* Operations */
        template<class OtherDerived>
        __host__ __device__ void assignTo(device_obj<LValueMatrix<OtherDerived>>& target) const;
        [[nodiscard]] __host__ __device__ inline RowVector row(size_t r);
        [[nodiscard]] __host__ __device__ inline const RowVector row(size_t r) const;
        [[nodiscard]] __host__ __device__ inline ColVector col(size_t c);
        [[nodiscard]] __host__ __device__ inline const ColVector col(size_t c) const;
        [[nodiscard]] __host__ __device__ inline BlockType rows(size_t fromRow, size_t rowCount);
        [[nodiscard]] __host__ __device__ inline const BlockType rows(size_t fromRow, size_t rowCount) const;
        [[nodiscard]] __host__ __device__ inline BlockType topRows(size_t to);
        [[nodiscard]] __host__ __device__ inline const BlockType topRows(size_t to) const;
        [[nodiscard]] __host__ __device__ inline BlockType bottomRows(size_t from);
        [[nodiscard]] __host__ __device__ inline const BlockType bottomRows(size_t from) const;
        [[nodiscard]] __host__ __device__ inline BlockType cols(size_t fromCol, size_t colCount);
        [[nodiscard]] __host__ __device__ inline const BlockType cols(size_t fromCol, size_t colCount) const;
        [[nodiscard]] __host__ __device__ inline BlockType leftCols(size_t to);
        [[nodiscard]] __host__ __device__ inline const BlockType leftCols(size_t to) const;
        [[nodiscard]] __host__ __device__ inline BlockType rightCols(size_t from);
        [[nodiscard]] __host__ __device__ inline const BlockType rightCols(size_t from) const;
        [[nodiscard]] __host__ __device__ inline BlockType topLeftCorner(size_t toRow, size_t toCol);
        [[nodiscard]] __host__ __device__ inline const BlockType topLeftCorner(size_t toRow, size_t toCol) const;
        [[nodiscard]] __host__ __device__ inline BlockType topLeftCorner(size_t to);
        [[nodiscard]] __host__ __device__ inline const BlockType topLeftCorner(size_t to) const;
        [[nodiscard]] __host__ __device__ inline BlockType topRightCorner(size_t toRow, size_t fromCol);
        [[nodiscard]] __host__ __device__ inline const BlockType topRightCorner(size_t toRow, size_t fromCol) const;
        [[nodiscard]] __host__ __device__ inline BlockType bottomLeftCorner(size_t fromRow, size_t toCol);
        [[nodiscard]] __host__ __device__ inline const BlockType bottomLeftCorner(size_t fromRow, size_t toCol) const;
        [[nodiscard]] __host__ __device__ inline BlockType bottomRightCorner(size_t fromRow, size_t fromCol);
        [[nodiscard]] __host__ __device__ inline const BlockType bottomRightCorner(size_t fromRow, size_t fromCol) const;
        [[nodiscard]] __host__ __device__ inline BlockType bottomRightCorner(size_t from);
        [[nodiscard]] __host__ __device__ inline const BlockType bottomRightCorner(size_t from) const;
        [[nodiscard]] __host__ __device__ inline BlockType block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount);
        [[nodiscard]] __host__ __device__ inline const BlockType block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const;
        /* Operations */
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return Base::getDerived().calc(row, col); }
        [[nodiscard]] __device__ ScalarType calcFromMajorMinor(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ device_obj<Transpose<Derived>> transpose() const noexcept;
        [[nodiscard]] __host__ __device__ device_obj<RValueFlatten<Derived>> flatten() const noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getColumn() const noexcept { return Base::getDerived().getColumn(); }
        [[nodiscard]] __host__ __device__ size_t getMaxMajor() const noexcept { return MatrixOption::getMaxMajor<device_obj<Derived>>(Base::getDerived()); }
        [[nodiscard]] __host__ __device__ size_t getMaxMinor() const noexcept { return MatrixOption::getMaxMinor<device_obj<Derived>>(Base::getDerived()); }
    };
}

#ifdef __CUDA_ARCH__
    #include "RValueMatrixImpl.cuh"
#endif
