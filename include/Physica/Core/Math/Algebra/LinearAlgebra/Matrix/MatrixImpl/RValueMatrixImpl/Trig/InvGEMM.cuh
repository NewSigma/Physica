/*
 * Copyright 2026 Weibo He.
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
    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
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
        void assign_mkl(Matrix auto& target) const noexcept;
        /* Getters */
        [[nodiscard]] size_t getRow() const noexcept { return getLHS().getRow(); }
        [[nodiscard]] size_t getCol() const noexcept { return getRHS().getCol(); }
        [[nodiscard]] auto&& getLHS(this auto&&) noexcept;
        [[nodiscard]] auto&& getRHS(this auto&&) noexcept;
    };

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    device_obj<GEMM<M1, M2>>::device_obj(Ref1 lhs, Ref2 rhs) : lhs(asStruct(lhs)), rhs(asStruct(rhs)) {}

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    void device_obj<GEMM<M1, M2>>::assign(Matrix auto& target) const {
        using Tm = decltype(std::declval<T>().toCUDA());
        if constexpr (Internal::isInvTrig<M1>())
            getRHS().assign(target);
        else
            getLHS().assign(target);

        auto& ctx = CUDAContext::getInstance();
        ctx.setPointerMode(false);
        constexpr auto Side = Internal::isInvTrig<M1>() ? CUBLAS_SIDE_LEFT : CUBLAS_SIDE_RIGHT;
        constexpr auto Uplo = []() consteval noexcept {
            if constexpr (Internal::isInvTrig<M1>())
                return Traits<M1>::Upper ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER;
            else
                return Traits<M2>::Upper ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER;
        }();
        constexpr auto Trans = CUBLAS_OP_N;
        constexpr auto Diag = []() consteval noexcept {
            if constexpr (Internal::isInvTrig<M1>())
                return Traits<M1>::Unit ? CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT;
            else
                return Traits<M2>::Unit ? CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT;
        }();
        const size_t m = getRow();
        const size_t n = getCol();
        const auto alpha = Tm(T(1));
        const auto* A = reinterpret_cast<const Tm*>([&]() noexcept {
            if constexpr (Internal::isInvTrig<M1>())
                return getLHS().getExpr().getExpr().data();
            else
                return getRHS().getExpr().getExpr().data();
        }());
        const size_t lda = Side == CUBLAS_SIDE_LEFT ? m : n;
        auto* B = reinterpret_cast<Tm*>(target.data());
        const size_t ldb = MatrixMajor::isColMatrix<M1>() ? m : n;
        if constexpr (Base::isComplex()) {
            if constexpr (T::Prec == Float32)
                check(cublasCtrsm_64(ctx, Side, Uplo, Trans, Diag, m, n, &alpha, A, lda, B, ldb));
            else
                check(cublasZtrsm_64(ctx, Side, Uplo, Trans, Diag, m, n, &alpha, A, lda, B, ldb));
        }
        else {
            if constexpr (T::Prec == Float32)
                check(cublasStrsm_64(ctx, Side, Uplo, Trans, Diag, m, n, &alpha, A, lda, B, ldb));
            else
                check(cublasDtrsm_64(ctx, Side, Uplo, Trans, Diag, m, n, &alpha, A, lda, B, ldb));
        }
    }

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    auto&& device_obj<GEMM<M1, M2>>::getLHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref1>(self.lhs.getDerived());
    }

    template<Matrix M1, Matrix M2> requires(Internal::isInvTrig<M1>() != Internal::isInvTrig<M2>())
    auto&& device_obj<GEMM<M1, M2>>::getRHS(this auto&& self) noexcept {
        return propagate_rvalue_reference<decltype(self), Ref2>(self.rhs.getDerived());
    }
}
