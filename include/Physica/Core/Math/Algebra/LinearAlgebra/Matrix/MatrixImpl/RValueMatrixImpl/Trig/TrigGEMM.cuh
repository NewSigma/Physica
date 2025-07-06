/*
 * Copyright 2025 Weibo He.
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

#include "Trig.cuh"

namespace Physica {
    template<Matrix M1, Matrix M2>
    class device_obj<GEMM<TrigUpper<M1>, M2>> : public device_obj<RValueMatrix<GEMM<TrigUpper<M1>, M2>>> {
        using host_obj = GEMM<TrigUpper<M1>, M2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
    public:
        using Base::isComplex;
    protected:
        using typename Base::T;
        using Tc = T::ComplexType;
        using Tm = std::conditional<isComplex, typename Tc::MKL_Complex, typename T::MachineType>::type;
    private:
        const device_obj<TrigUpper<M1>>& mat1;
        const device_obj<M2>& mat2;
    public:
        device_obj(const device_obj<TrigUpper<M1>>& mat1_, const device_obj<M2>& mat2_);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return mat1.getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return mat2.getCol(); }
        [[nodiscard]] __host__ __device__ const auto& getLHS() const noexcept { return mat1; }
        [[nodiscard]] __host__ __device__ const auto& getRHS() const noexcept { return mat2; }
    };

    template<Matrix M1, Matrix M2>
    device_obj<GEMM<TrigUpper<M1>, M2>>::device_obj(const device_obj<TrigUpper<M1>>& mat1_, const device_obj<M2>& mat2_) : mat1(mat1_), mat2(mat2_) {}

    template<Matrix M1, Matrix M2>
    void device_obj<GEMM<TrigUpper<M1>, M2>>::assign(Matrix auto& target) const {
        using M = std::remove_cvref_t<decltype(target)>;
        constexpr auto Layout = MatrixOption::isRowMatrix<M>() ? CblasRowMajor : CblasColMajor;
        constexpr auto Side = CUBLAS_SIDE_LEFT;
        constexpr auto Uplo = CUBLAS_FILL_MODE_UPPER;
        constexpr auto TransA = CUBLAS_OP_N;
        constexpr auto Diag = CUBLAS_DIAG_NON_UNIT;
        const M buffer = mat1;
        target = mat2;

        const size_t m = getRow();
        const size_t n = getCol();
        const Tm alpha = 1;
        const auto* A = reinterpret_cast<const Tm*>(buffer.data());
        const size_t lda = Side == CUBLAS_SIDE_LEFT ? m : n;
        auto* B = reinterpret_cast<Tm*>(target.data());
        const size_t ldb = Layout == CblasColMajor ? m : n;

        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        if constexpr (isComplex) {
            if constexpr (T::Prec == Float32)
                check(cublasCtrmm_64(ctx, Side, Uplo, TransA, Diag, m, n, &alpha, A, lda, B, ldb, B, ldb));
            else
                check(cublasZtrmm_64(ctx, Side, Uplo, TransA, Diag, m, n, &alpha, A, lda, B, ldb, B, ldb));
        }
        else {
            if constexpr (T::Prec == Float32)
                check(cublasStrmm_64(ctx, Side, Uplo, TransA, Diag, m, n, &alpha, A, lda, B, ldb, B, ldb));
            else
                check(cublasDtrmm_64(ctx, Side, Uplo, TransA, Diag, m, n, &alpha, A, lda, B, ldb, B, ldb));
        }
    }
}
