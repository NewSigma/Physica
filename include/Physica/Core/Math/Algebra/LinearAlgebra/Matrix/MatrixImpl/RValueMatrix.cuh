/*
 * Copyright 2022-2024 Weibo He.
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

#include "RValueMatrixImpl/RMatrixBlock.cuh"

namespace Physica::Core {
    template<class Derived>
    class device_obj<RValueMatrix<Derived>> : public CRTPBase<device_obj<RValueMatrix<Derived>>> {
        static_assert(!is_device_obj<Derived>::value, "[Error]: Nested device_obj is unnecessary");
        using host_obj = RValueMatrix<Derived>;
        using This = device_obj<host_obj>;
        using Base = CRTPBase<This>;
        using TraitsType = Traits<device_obj<Derived>>;
        using RowVector = device_obj<RMatrixBlock<Derived, 1, Dynamic>>;
        using ColVector = device_obj<RMatrixBlock<Derived, Dynamic, 1>>;
        using BlockType = device_obj<RMatrixBlock<Derived>>;
    public:
        using ScalarType = TraitsType::ScalarType;
        using ValueType = ScalarType::ValueType;
        constexpr static int Option = TraitsType::Option;
        constexpr static size_t RowAtCompile = TraitsType::RowAtCompile;
        constexpr static size_t ColAtCompile = TraitsType::ColAtCompile;
        constexpr static size_t SizeAtCompile = TraitsType::SizeAtCompile;
        constexpr static bool isReverseDiff = ScalarType::isReverseDiff;
        constexpr static bool isComplex = ScalarType::isComplex;
        constexpr static size_t MaxThreadPerBlock = 256;
    public:
        ~device_obj() = default;
        /* Operations */
        template<Matrix M>
        __host__ __device__ void assignTo(device_obj<LValueMatrix<M>>& target) const;
        [[nodiscard]] __host__ __device__ inline auto row(size_t r) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto row(size_t r) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto col(size_t c) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto col(size_t c) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto rows(size_t fromRow, size_t rowCount) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto rows(size_t fromRow, size_t rowCount) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto topRows(size_t to) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto topRows(size_t to) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto bottomRows(size_t from) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto bottomRows(size_t from) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto cols(size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto cols(size_t fromCol, size_t colCount) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto leftCols(size_t to) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto leftCols(size_t to) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto rightCols(size_t from) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto rightCols(size_t from) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto topLeftCorner(size_t toRow, size_t toCol) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto topLeftCorner(size_t toRow, size_t toCol) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto topLeftCorner(size_t to) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto topLeftCorner(size_t to) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto topRightCorner(size_t toRow, size_t fromCol) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto topRightCorner(size_t toRow, size_t fromCol) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto bottomLeftCorner(size_t fromRow, size_t toCol) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto bottomRightCorner(size_t fromRow, size_t fromCol) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto bottomRightCorner(size_t from) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto bottomRightCorner(size_t from) const noexcept;
        [[nodiscard]] __host__ __device__ inline auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] __host__ __device__ inline const auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept;
        /* Operations */
        template<Side Owner = GetSide()>
        [[nodiscard]] __device__ ScalarType calc(size_t row, size_t col) const { return Base::getDerived().template calc<Owner>(row, col); }
        [[nodiscard]] __device__ inline ScalarType calcFromMajorMinor(size_t row, size_t col) const;
        [[nodiscard]] __host__ __device__ auto transpose() const noexcept;
        [[nodiscard]] __host__ __device__ auto flatten() const noexcept;
        /* Getters */
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getRow() const noexcept { return Base::getDerived().template getRow<Owner>(); }
        template<Side Owner = GetSide()>
        [[nodiscard]] __host__ __device__ size_t getCol() const noexcept { return Base::getDerived().template getCol<Owner>(); }
        [[nodiscard]] __host__ __device__ size_t getMaxMajor() const noexcept { return MatrixOption::getMaxMajor<device_obj<Derived>>(Base::getDerived()); }
        [[nodiscard]] __host__ __device__ size_t getMaxMinor() const noexcept { return MatrixOption::getMaxMinor<device_obj<Derived>>(Base::getDerived()); }
        /* Static members */
        [[nodiscard]] __host__ __device__ static size_t rowFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::rowFromMajorMinor<device_obj<Derived>>(major, minor); }
        [[nodiscard]] __host__ __device__ static size_t colFromMajorMinor(size_t major, size_t minor) noexcept { return MatrixOption::colFromMajorMinor<device_obj<Derived>>(major, minor); }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        /* Operators */
        This& operator=(const This&) = default;
        This& operator=(This&&) noexcept = default;
    };
}

namespace Physica {
    template<class T>
    class Traits<device_obj<RValueMatrix<T>>> {
    public:
        using Derived = device_obj<T>;
    };
}

#include "RValueMatrixImpl/RValueMatrixImpl.cuh"
#include "Transpose.cuh"
#include "MatrixExpr.cuh"
#include "MatrixProduct/GEMM.cuh"
#include "MatrixProduct/GEMV.cuh"
#include "MatrixProduct/GEVM.cuh"
