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

#include "RValueMatrix.cuh"
#include "LValueMatrixImpl/LMatrixBlock.cuh"

namespace Physica {
    template<class Derived>
    class device_obj<LValueMatrix<Derived>> : public device_obj<RValueMatrix<Derived>> {
        using This = device_obj<LValueMatrix<Derived>>;
        using Base = device_obj<RValueMatrix<Derived>>;
        using RowVector = device_obj<LMatrixBlock<Derived, 1, Dynamic>>;
        using ColVector = device_obj<LMatrixBlock<Derived, Dynamic, 1>>;
        using BlockType = device_obj<LMatrixBlock<Derived>>;
    public:
        using typename Base::ScalarType;
        using Base::RowAtCompile;
        using Base::ColAtCompile;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
        using PtrTy = T::PtrTy;
        using ConstPtrTy = T::ConstPtrTy;
        using RefTy = T::RefTy;
        using ConstRefTy = T::ConstRefTy;
    public:
        /* Operators */
        This& operator=(const This& m) = delete;
        This& operator=(This&& m) noexcept = delete;

        template<Scalar U>
        __host__ __device__ device_obj<Derived>& operator=(const U& x);
        __host__ __device__ void operator+=(const Scalar auto& x);
        __host__ __device__ void operator-=(const Scalar auto& x);
        __host__ __device__ void operator*=(const Scalar auto& x);
        __host__ __device__ void operator/=(const Scalar auto& x);

        __host__ __device__ void operator+=(const Vector auto& v) { Base::getDerived() = Base::getDerived() + v; }
        __host__ __device__ void operator-=(const Vector auto& v) { Base::getDerived() = Base::getDerived() - v; }

        __host__ __device__ device_obj<Derived>& operator=(const Matrix auto& m);
        __host__ __device__ void operator+=(const Matrix auto& m);
        __host__ __device__ void operator-=(const Matrix auto& m);

        [[nodiscard]] __device__ RefTy operator()(size_t row, size_t col);
        [[nodiscard]] __device__ ConstRefTy operator()(size_t row, size_t col) const;
        /* Operations */
        [[nodiscard]] __device__ ConstRefTy calc(size_t row, size_t col) const { return operator()(row, col); }
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }

        void reverse(const Matrix auto& grad) const noexcept;

        [[nodiscard]] __host__ __device__ auto row(size_t r) noexcept;
        [[nodiscard]] __host__ __device__ const auto row(size_t r) const noexcept;
        [[nodiscard]] __host__ __device__ auto col(size_t c) noexcept;
        [[nodiscard]] __host__ __device__ const auto col(size_t c) const noexcept;
        [[nodiscard]] __host__ __device__ auto rows(size_t fromRow, size_t rowCount) noexcept;
        [[nodiscard]] __host__ __device__ const auto rows(size_t fromRow, size_t rowCount) const noexcept;
        [[nodiscard]] __host__ __device__ auto topRows(size_t to) noexcept;
        [[nodiscard]] __host__ __device__ const auto topRows(size_t to) const noexcept;
        [[nodiscard]] __host__ __device__ auto bottomRows(size_t from) noexcept;
        [[nodiscard]] __host__ __device__ const auto bottomRows(size_t from) const noexcept;
        [[nodiscard]] __host__ __device__ auto cols(size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] __host__ __device__ const auto cols(size_t fromCol, size_t colCount) const noexcept;
        [[nodiscard]] __host__ __device__ auto leftCols(size_t to) noexcept;
        [[nodiscard]] __host__ __device__ const auto leftCols(size_t to) const noexcept;
        [[nodiscard]] __host__ __device__ auto rightCols(size_t from) noexcept;
        [[nodiscard]] __host__ __device__ const auto rightCols(size_t from) const noexcept;
        [[nodiscard]] __host__ __device__ auto topLeftCorner(size_t toRow, size_t toCol) noexcept;
        [[nodiscard]] __host__ __device__ const auto topLeftCorner(size_t toRow, size_t toCol) const noexcept;
        [[nodiscard]] __host__ __device__ auto topLeftCorner(size_t to) noexcept;
        [[nodiscard]] __host__ __device__ const auto topLeftCorner(size_t to) const noexcept;
        [[nodiscard]] __host__ __device__ auto topRightCorner(size_t toRow, size_t fromCol) noexcept;
        [[nodiscard]] __host__ __device__ const auto topRightCorner(size_t toRow, size_t fromCol) const noexcept;
        [[nodiscard]] __host__ __device__ auto bottomLeftCorner(size_t fromRow, size_t toCol) noexcept;
        [[nodiscard]] __host__ __device__ const auto bottomLeftCorner(size_t fromRow, size_t toCol) const noexcept;
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(size_t fromRow, size_t fromCol) noexcept;
        [[nodiscard]] __host__ __device__ const auto bottomRightCorner(size_t fromRow, size_t fromCol) const noexcept;
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(size_t from) noexcept;
        [[nodiscard]] __host__ __device__ const auto bottomRightCorner(size_t from) const noexcept;
        [[nodiscard]] __host__ __device__ auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;
        [[nodiscard]] __host__ __device__ const auto block(size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) const noexcept;

        [[nodiscard]] __host__ __device__ auto flatten();
        [[nodiscard]] __host__ __device__ const auto flatten() const;

        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();
        template<RNG R>
        void random_any(auto& distribution);
        /* Getters */
        using Base::getRow;
        using Base::getCol;
        [[nodiscard]] __host__ __device__ PtrTy data_ptr(size_t row, size_t col) noexcept;
        [[nodiscard]] __host__ __device__ ConstPtrTy data_ptr(size_t row, size_t col) const noexcept;
        [[nodiscard]] __device__ RefTy refFromMajorMinor(size_t major, size_t minor) noexcept;
        [[nodiscard]] __device__ ConstRefTy refFromMajorMinor(size_t major, size_t minor) const noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "LValueMatrixImpl/LValueMatrixImpl.cuh"
#include "LValueMatrixImpl/Flatten.cuh"
#include "LValueMatrixImpl/ReshapedVector.cuh"
