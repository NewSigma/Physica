/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "CompactMatrix.h"
#include "StridedMatrix.cuh"
#include "CompactMatrixImpl/CompactMatrixBlock.cuh"

namespace Physica {
    template<class Derived>
    class device_obj<CompactMatrix<Derived>> : public device_obj<StridedMatrix<Derived>> {
        using host_obj = CompactMatrix<Derived>;
        using This = device_obj<host_obj>;
        using Base = device_obj<StridedMatrix<Derived>>;
        template<class Self, size_t Row, size_t Col>
        using Block = device_obj<typename CompactMatrix<Derived>::template Block<Self, Row, Col>>;
    public:
        using typename Base::ScalarType;
    protected:
        using typename Base::T;
    private:
        constexpr static bool isRowMatrix = MatrixMajor::isRowMatrix<Derived>();
        constexpr static bool isColMatrix = MatrixMajor::isColMatrix<Derived>();
    public:
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Operations */
        template<Matrix M> void toHost(CompactMatrix<M>& obj) const;
        template<Matrix M> void toHostAsync(CompactMatrix<M>& obj) const;

        [[nodiscard]] __host__ __device__ auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] __host__ __device__ auto col(this auto&&, size_t c) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ auto rows(this auto&&, size_t fromRow, size_t rowCount) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ auto topRows(this auto&&, size_t to) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] __host__ __device__ auto bottomRows(this auto&&, size_t from) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto cols(this auto&&, size_t fromCol, size_t colCount) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto leftCols(this auto&&, size_t to) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto rightCols(this auto&&, size_t from) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto topLeftCorner(this auto&&, size_t toRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto topLeftCorner(this auto&&, size_t to) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto topRightCorner(this auto&&, size_t toRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto bottomLeftCorner(this auto&&, size_t fromRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(this auto&&, size_t fromRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto bottomRightCorner(this auto&&, size_t from) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] __host__ __device__ auto block(this auto&&, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;

        void zeros();
        template<RNG R>
        void random_uniform();
        template<RNG R>
        void random_normal();

        void read(const auto& obj);
        /* Getters */
        [[nodiscard]] __host__ __device__ auto data(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getRowStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getColStride() const noexcept;
    protected:
        device_obj() = default;
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
    };
}

#include "CompactMatrixImpl/CompactMatrixImpl.cuh"
#include "CompactMatrixImpl/Flatten.cuh"
#include "CompactMatrixImpl/Transpose.cuh"
