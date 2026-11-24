/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
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
        using Base::isReverseDiff;
    protected:
        using typename Base::T;
        using typename Base::Tr;
        using typename Base::Tv;
        using typename Base::Trv;
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
        [[nodiscard]] __device__ T calc(size_t row, size_t col) const;
        [[nodiscard]] __device__ T calc(size_t row, size_t col, instanceof_x<ThreadBlock> auto) const;

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

        __host__ __device__ void zeros() noexcept;
        __host__ __device__ void zero_grad() noexcept;
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
#include "LValueMatrixImpl/MainDiag.cuh"
#include "LValueMatrixImpl/OffsetDiag.cuh"
#include "LValueMatrixImpl/Flatten.cuh"
#include "LValueMatrixImpl/ReshapedVector.cuh"
