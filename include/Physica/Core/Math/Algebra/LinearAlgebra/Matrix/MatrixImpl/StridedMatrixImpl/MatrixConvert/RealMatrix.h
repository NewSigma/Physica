/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/StridedMatrix.h"

namespace Physica {
    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    class RealMatrix<M> : public StridedMatrix<RealMatrix<M>> {
        using This = RealMatrix<M>;
        using Base = StridedMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit RealMatrix(M&& mat_);
        RealMatrix(const This&) = default;
        RealMatrix(This&&) noexcept = default;
        ~RealMatrix() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        using Base::operator=;
        /* Operations */
        using Base::resize;
        void resize(size_t row, size_t col);
        /* Getters */
        [[nodiscard]] auto data_handle(this auto&& self) noexcept;
        [[nodiscard]] constexpr size_t getRow() const noexcept;
        [[nodiscard]] constexpr size_t getCol() const noexcept;
        [[nodiscard]] constexpr size_t getOrder() const noexcept;
        [[nodiscard]] constexpr size_t getRowStride() const noexcept;
        [[nodiscard]] constexpr size_t getColStride() const noexcept;
        /* Static members */
        [[nodiscard]] __host__ __device__ consteval static bool isStaticSquare() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getRowAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static int getMajor() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getRowStrideAtCompile() noexcept;
        [[nodiscard]] __host__ __device__ consteval static size_t getColStrideAtCompile() noexcept;
    };

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    RealMatrix<M>::RealMatrix(M&& mat_) : mat(std::forward<M>(mat_)) {}

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    void RealMatrix<M>::resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    auto RealMatrix<M>::data_handle(this auto&& self) noexcept {
        return self.mat.data_handle()->real_ptr();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t RealMatrix<M>::getRow() const noexcept {
        return mat.getRow();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t RealMatrix<M>::getCol() const noexcept {
        return mat.getCol();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t RealMatrix<M>::getOrder() const noexcept {
        return mat.getOrder();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t RealMatrix<M>::getRowStride() const noexcept {
        return 2 * mat.getRowStride();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t RealMatrix<M>::getColStride() const noexcept {
        return 2 * mat.getColStride();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval bool RealMatrix<M>::isStaticSquare() noexcept {
        return std::remove_cvref_t<M>::isStaticSquare();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t RealMatrix<M>::getRowAtCompile() noexcept {
        return std::remove_cvref_t<M>::getRowAtCompile();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t RealMatrix<M>::getColAtCompile() noexcept {
        return std::remove_cvref_t<M>::getColAtCompile();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval int RealMatrix<M>::getMajor() noexcept {
        return std::remove_cvref_t<M>::getMajor();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t RealMatrix<M>::getRowStrideAtCompile() noexcept {
        return 2 * std::remove_cvref_t<M>::getRowStrideAtCompile();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t RealMatrix<M>::getColStrideAtCompile() noexcept {
        return 2 * std::remove_cvref_t<M>::getColStrideAtCompile();
    }
}
