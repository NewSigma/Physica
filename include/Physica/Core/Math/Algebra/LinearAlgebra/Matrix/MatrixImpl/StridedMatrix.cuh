/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "LValueMatrix.cuh"
#include "StridedMatrix.h"

namespace Physica {
    template<class Derived>
    class device_obj<StridedMatrix<Derived>> : public device_obj<LValueMatrix<Derived>> {
        using host_obj = StridedMatrix<Derived>;
        using Base = device_obj<LValueMatrix<Derived>>;
        using This = device_obj<host_obj>;
    protected:
        using typename Base::T;
    public:
        constexpr ~device_obj() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;
        /* Operations */
        __host__ __device__ void assign(this const auto& self, Matrix auto&& target);
        __device__ void assign(this const auto&, Matrix auto&& target, instanceof_x<ThreadBlock> auto block);

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

        __host__ __device__ void zeros() noexcept;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getRowStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getColStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getMajorStride() const noexcept;
        [[nodiscard]] __host__ __device__ constexpr size_t getMinorStride() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() noexcept;
        [[nodiscard]] __host__ __device__ auto data_handle() const noexcept;
        [[nodiscard]] __host__ __device__ auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static size_t getRowStrideAtCompile() noexcept { return host_obj::getRowStrideAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept { return host_obj::getColStrideAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getMajorStrideAtCompile() noexcept { return host_obj::getMajorStrideAtCompile(); }
        [[nodiscard]] __host__ __device__ consteval static size_t getMinorStrideAtCompile() noexcept { return host_obj::getMinorStrideAtCompile(); }
    protected:
        constexpr device_obj() = default;
        constexpr device_obj(const This&) = default;
        constexpr device_obj(This&&) noexcept = default;
    };
}

#include "StridedMatrixImpl/StridedMatrixImpl.cuh"
#include "StridedMatrixImpl/StridedMatrixBlock.cuh"
#include "StridedMatrixImpl/MainDiag.cuh"
#include "StridedMatrixImpl/OffsetDiag.cuh"
#include "StridedMatrixImpl/MatrixConvert/RealMatrix.cuh"
#include "StridedMatrixImpl/MatrixConvert/ImagMatrix.cuh"
