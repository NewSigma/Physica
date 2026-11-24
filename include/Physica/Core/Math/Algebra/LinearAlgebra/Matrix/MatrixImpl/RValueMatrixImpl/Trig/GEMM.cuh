/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "MatrixTrig.cuh"

namespace Physica {
    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
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
        PlainStruct<add_device_obj_t<std::remove_reference_t<M1>>> lhs;
        PlainStruct<add_device_obj_t<std::remove_reference_t<M2>>> rhs;
    public:
        device_obj(Ref1 lhs, Ref2 rhs);
        device_obj(const This&) = default;
        device_obj(This&&) noexcept = default;
        ~device_obj() = default;
        /* Operators */
        This& operator=(const This&) = delete;
        This& operator=(This&&) noexcept = delete;
        /* Operations */
        void assign(Matrix auto& target) const;
        /* Getters */
        [[nodiscard]] __host__ __device__ constexpr size_t getRow() const { return getLHS().getRow(); }
        [[nodiscard]] __host__ __device__ constexpr size_t getCol() const { return getRHS().getCol(); }
        [[nodiscard]] __host__ __device__ auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] __host__ __device__ auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    device_obj<GEMM<M1, M2>>::device_obj(Ref1 lhs, Ref2 rhs) : lhs(asStruct(lhs)), rhs(asStruct(rhs)) {}

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    void device_obj<GEMM<M1, M2>>::assign(Matrix auto& target) const {
        using M = std::remove_cvref_t<decltype(target)>;
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr bool TrigLHS = instanceof_tx<M1, MatrixTrig>;
        using Trig = std::conditional_t<TrigLHS, M1, M2>;
        using Data = std::remove_cvref_t<decltype(std::declval<const Trig&>().getExpr())>;
        static_assert(Data::getMinorStrideAtCompile() == 1 && M::getMinorStrideAtCompile() == 1,
                      "[Error]: cuBLAS expects a contiguous minor direction");
        constexpr auto Side = TrigLHS ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT;
        constexpr auto TransA = MatrixMajor::isSameMajor<Trig, M>() ? CUBLAS_OP_N : CUBLAS_OP_T;
        constexpr auto UploNoTrans = Traits<Trig>::Upper ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER;
        constexpr auto Uplo = []() consteval noexcept {
            if constexpr (TransA == CUBLAS_OP_T)
                return UploNoTrans == CUBLAS_FILL_MODE_UPPER ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
            return UploNoTrans;
        }();
        constexpr auto Diag = Traits<Trig>::Unit ? CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT;

        [this]() -> auto& {
            if constexpr (TrigLHS)
                return getRHS();
            else
                return getLHS();
        }().assign(target);

        const size_t m = getRow();
        const size_t n = getCol();
        const Tm alpha = 1;
        const auto* A = reinterpret_cast<const Tm*>([this]() noexcept {
            if constexpr (TrigLHS)
                return getLHS().getExpr().data_handle();
            else
                return getRHS().getExpr().data_handle();
        }());
        const size_t lda = [this]() noexcept {
            if constexpr (TrigLHS)
                return getLHS().getExpr().getMajorStride();
            else
                return getRHS().getExpr().getMajorStride();
        }();
        auto* B = reinterpret_cast<Tm*>(target.data_handle());
        const size_t ldb = target.getMajorStride();

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

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    __host__ __device__ auto&& device_obj<GEMM<M1, M2>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref1>(self.lhs.getDerived());
    }

    template<Matrix M1, Matrix M2> requires(instanceof_tx<M1, MatrixTrig> || instanceof_tx<M2, MatrixTrig>)
    __host__ __device__ auto&& device_obj<GEMM<M1, M2>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref2>(self.rhs.getDerived());
    }
}
