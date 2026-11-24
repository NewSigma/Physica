/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixImpl/StridedMatrix.h"

namespace Physica {
    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    class ImagMatrix<M> : public StridedMatrix<ImagMatrix<M>> {
        using This = ImagMatrix<M>;
        using Base = StridedMatrix<This>;
    protected:
        using typename Base::T;
    private:
        decay_rvalue_t<M> mat;
    public:
        explicit ImagMatrix(M&& mat_);
        ImagMatrix(const This&) = default;
        ImagMatrix(This&&) noexcept = default;
        ~ImagMatrix() = default;
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
    ImagMatrix<M>::ImagMatrix(M&& mat_) : mat(std::forward<M>(mat_)) {}

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    void ImagMatrix<M>::resize([[maybe_unused]] size_t row, [[maybe_unused]] size_t col) {
        assert(row == getRow());
        assert(col == getCol());
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    auto ImagMatrix<M>::data_handle(this auto&& self) noexcept {
        return self.mat.data_handle()->imag_ptr();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t ImagMatrix<M>::getRow() const noexcept {
        return mat.getRow();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t ImagMatrix<M>::getCol() const noexcept {
        return mat.getCol();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t ImagMatrix<M>::getOrder() const noexcept {
        return mat.getOrder();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t ImagMatrix<M>::getRowStride() const noexcept {
        return 2 * mat.getRowStride();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    constexpr size_t ImagMatrix<M>::getColStride() const noexcept {
        return 2 * mat.getColStride();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval bool ImagMatrix<M>::isStaticSquare() noexcept {
        return std::remove_cvref_t<M>::isStaticSquare();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t ImagMatrix<M>::getRowAtCompile() noexcept {
        return std::remove_cvref_t<M>::getRowAtCompile();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t ImagMatrix<M>::getColAtCompile() noexcept {
        return std::remove_cvref_t<M>::getColAtCompile();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval int ImagMatrix<M>::getMajor() noexcept {
        return std::remove_cvref_t<M>::getMajor();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t ImagMatrix<M>::getRowStrideAtCompile() noexcept {
        return 2 * std::remove_cvref_t<M>::getRowStrideAtCompile();
    }

    template<class M> requires(std::remove_cvref_t<M>::isStrided() && std::remove_cvref_t<M>::isComplex())
    __host__ __device__ consteval size_t ImagMatrix<M>::getColStrideAtCompile() noexcept {
        return 2 * std::remove_cvref_t<M>::getColStrideAtCompile();
    }
}
