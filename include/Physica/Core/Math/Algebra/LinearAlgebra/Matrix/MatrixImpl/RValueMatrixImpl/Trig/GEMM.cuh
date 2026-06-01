/*
 * Copyright 2025-2026 Weibo He.
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

#include "MatrixTrig.cuh"

namespace Physica {
    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    class device_obj<GEMM<M1, M2>> : public device_obj<RValueMatrix<GEMM<M1, M2>>> {
        using host_obj = GEMM<M1, M2>;
        using This = device_obj<host_obj>;
        using Base = device_obj<RValueMatrix<host_obj>>;
        using Ref1 = add_device_obj<M1>::type;
        using Ref2 = add_device_obj<M2>::type;
    protected:
        using typename Base::T;
        using typename Base::Tc;
    private:
        PlainStruct<add_device_obj_t<std::remove_reference_t<M1>>> trig;
        PlainStruct<add_device_obj_t<std::remove_reference_t<M2>>> rhs;
    public:
        device_obj(Ref1 trig, Ref2 rhs);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ size_t getRow() const { return getLHS().getRow(); }
        [[nodiscard]] __host__ __device__ size_t getCol() const { return getRHS().getCol(); }
        [[nodiscard]] __host__ __device__ auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    device_obj<GEMM<M1, M2>>::device_obj(Ref1 trig, Ref2 rhs) : trig(asStruct(trig)), rhs(asStruct(rhs)) {}

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    void device_obj<GEMM<M1, M2>>::assign(Matrix auto& target) const {
        using M = std::remove_cvref_t<decltype(target)>;
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr auto Side = CUBLAS_SIDE_LEFT;
        constexpr auto Uplo = Traits<M1>::Upper ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER;
        constexpr auto TransA = CUBLAS_OP_N;
        constexpr auto Diag = Traits<M1>::Unit ? CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT;
        const M buffer = getLHS();
        getRHS().assign(target);

        const size_t m = getRow();
        const size_t n = getCol();
        const Tm alpha = 1;
        const auto* A = reinterpret_cast<const Tm*>(buffer.data());
        const size_t lda = Side == CUBLAS_SIDE_LEFT ? m : n;
        auto* B = reinterpret_cast<Tm*>(target.data());
        const size_t ldb = MatrixMajor::isColMatrix<M>() ? m : n;

        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        if constexpr (Base::isComplex()) {
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

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    __host__ __device__ auto&& device_obj<GEMM<M1, M2>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref1>(self.trig.getDerived());
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig>)
    __host__ __device__ auto&& device_obj<GEMM<M1, M2>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref2>(self.rhs.getDerived());
    }
}
