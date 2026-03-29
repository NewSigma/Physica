/*
 * Copyright 2022-2026 Weibo He.
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
    public:
        using typename Base::ScalarType;
        using Base::RowAtCompile;
        using Base::ColAtCompile;
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
    public:
        /* Operators */
        This& operator=(const This& m) = delete;
        This& operator=(This&& m) noexcept = delete;

        __host__ __device__ auto operator=(Scalar auto x) -> device_obj<Derived>&;
        __host__ __device__ void operator+=(Scalar auto x);
        __host__ __device__ void operator-=(Scalar auto x);
        __host__ __device__ void operator*=(Scalar auto x);
        __host__ __device__ void operator/=(Scalar auto x);

        __host__ __device__ void operator+=(const Vector auto& v) { Base::getDerived() = Base::getDerived() + v; }
        __host__ __device__ void operator-=(const Vector auto& v) { Base::getDerived() = Base::getDerived() - v; }

        __host__ __device__ device_obj<Derived>& operator=(const Matrix auto& m);
        __host__ __device__ void operator+=(const Matrix auto& m);
        __host__ __device__ void operator-=(const Matrix auto& m);

        [[nodiscard]] __device__ decltype(auto) operator[](this auto&& self, size_t row, size_t col);
        /* Operations */
        [[nodiscard]] __device__ decltype(auto) calc(size_t row, size_t col) const { return operator[](row, col); }
        [[nodiscard]] __device__ Tv calc_value(size_t row, size_t col) const { return calc(row, col).value(); }

        void reverse(const Matrix auto& grad) const noexcept;

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
        [[nodiscard]] __host__ __device__ auto diag(this auto&&) noexcept;

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
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        [[nodiscard]] __device__ decltype(auto) refFromMajorMinor(this auto&&, size_t major, size_t minor) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isLValueMatrix() noexcept { return true; }
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "LValueMatrixImpl/LValueMatrixImpl.cuh"
#include "LValueMatrixImpl/DiagVector.cuh"
#include "LValueMatrixImpl/Flatten.cuh"
#include "LValueMatrixImpl/ReshapedVector.cuh"
