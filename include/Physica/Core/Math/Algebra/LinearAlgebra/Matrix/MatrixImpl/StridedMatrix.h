/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "LValueMatrix.h"

namespace Physica {
    template<class Derived>
    class StridedMatrix : public LValueMatrix<Derived> {
        using Base = LValueMatrix<Derived>;
        using This = StridedMatrix<Derived>;
    protected:
        using typename Base::T;
    public:
        constexpr ~StridedMatrix() = default;
        /* Operators */
        This& operator=(const This& obj) = delete;
        This& operator=(This&& obj) noexcept = delete;
        using Base::operator=;

        [[nodiscard]] bool operator==(this const auto& self, const Matrix auto& other) noexcept;
        /* Operations */
        template<ExecutePolicy P = Sequential>
        void assign(Matrix auto&& target) const noexcept;

        [[nodiscard]] constexpr auto row(this auto&&, size_t r) noexcept;
        [[nodiscard]] constexpr auto col(this auto&&, size_t c) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] auto rows(this auto&&, size_t fromRow, size_t rowCount) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] auto topRows(this auto&&, size_t to) noexcept;
        template<size_t Row = Dynamic>
        [[nodiscard]] auto bottomRows(this auto&&, size_t from) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] auto cols(this auto&&, size_t fromCol, size_t colCount) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] auto leftCols(this auto&&, size_t to) noexcept;
        template<size_t Col = Dynamic>
        [[nodiscard]] auto rightCols(this auto&&, size_t from) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t toRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto topLeftCorner(this auto&&, size_t to) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto topRightCorner(this auto&&, size_t toRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto bottomLeftCorner(this auto&&, size_t fromRow, size_t toCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t fromRow, size_t fromCol) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto bottomRightCorner(this auto&&, size_t from) noexcept;
        template<size_t Row = Dynamic, size_t Col = Dynamic>
        [[nodiscard]] auto block(this auto&&, size_t fromRow, size_t rowCount, size_t fromCol, size_t colCount) noexcept;

        void zeros() noexcept;
        /* Getters */
        [[nodiscard]] constexpr size_t getRowStride() const noexcept;
        [[nodiscard]] constexpr size_t getColStride() const noexcept;
        [[nodiscard]] constexpr size_t getMajorStride() const noexcept;
        [[nodiscard]] constexpr size_t getMinorStride() const noexcept;
        [[nodiscard]] auto data_handle() noexcept;
        [[nodiscard]] auto data_handle() const noexcept;
        [[nodiscard]] constexpr auto data_ptr(this auto&&, size_t row, size_t col) noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStrided() noexcept { return true; }
        [[nodiscard]] __host__ __device__ consteval static size_t getRowStrideAtCompile() noexcept { return Dynamic; }
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept { return Dynamic; }
        [[nodiscard]] __host__ __device__ consteval static size_t getMajorStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getMinorStrideAtCompile() noexcept;
    protected:
        constexpr StridedMatrix();
        constexpr StridedMatrix(const This&) = default;
        constexpr StridedMatrix(This&&) noexcept = default;
    private:
        template<class M1, class M2>
        [[nodiscard]] consteval static bool isMinorCompact() noexcept;
        /* Friends */
        friend class device_obj<This>;
    };
}

#include "StridedMatrixImpl/StridedMatrixBlock.h"
#include "StridedMatrixImpl/StridedMatrixImpl.h"
#include "StridedMatrixImpl/MainDiag.h"
#include "StridedMatrixImpl/OffsetDiag.h"
#include "StridedMatrixImpl/MatrixConvert/RealMatrix.h"
#include "StridedMatrixImpl/MatrixConvert/ImagMatrix.h"
