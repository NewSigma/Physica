/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Core/Exception/MKL/Lapack.h"
#include "DenseQR.h"

namespace Physica {
    template<Scalar T, bool Pivot>
    void DenseQR<T, Pivot>::compute_mkl(const Matrix auto& source) {
        static_assert(T::Prec == Float32 || T::Prec == Float64);
        using Tm = decltype(std::declval<T>().toMKL());
        assert(getRow() == source.getRow());
        assert(getCol() == source.getCol());
        working = source;

        constexpr int Layout = LAPACK_COL_MAJOR;
        const size_t m = getRow();
        const size_t n = getCol();
        const size_t lda = m;
        auto* a = reinterpret_cast<Tm*>(working.data());
        auto* tau = reinterpret_cast<Tm*>(taus.data());

        if constexpr (Pivot) {
            auto* jpvt = reinterpret_cast<MKL_INT64*>(perm.getIndices().data());
            if constexpr (isComplex) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgeqp3_64(Layout, m, n, a, lda, jpvt, tau));
                else
                    check_lapack(LAPACKE_zgeqp3_64(Layout, m, n, a, lda, jpvt, tau));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgeqp3_64(Layout, m, n, a, lda, jpvt, tau));
                else
                    check_lapack(LAPACKE_dgeqp3_64(Layout, m, n, a, lda, jpvt, tau));
            }

            for (auto& index : perm.getIndices())
                index -= 1;
        }
        else {
            if constexpr (isComplex) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgeqrf_64(Layout, m, n, a, lda, tau));
                else
                    check_lapack(LAPACKE_zgeqrf_64(Layout, m, n, a, lda, tau));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgeqrf_64(Layout, m, n, a, lda, tau));
                else
                    check_lapack(LAPACKE_dgeqrf_64(Layout, m, n, a, lda, tau));
            }
        }
    }

    template<Scalar T, bool Pivot>
    auto DenseQR<T, Pivot>::getMatrixQ_mkl() const -> MatrixND<T> {
        static_assert(T::Prec == Float32 || T::Prec == Float64);
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr int Layout = LAPACK_COL_MAJOR;
        const size_t m = getRow();
        const size_t n = m;
        const size_t k = getCol();
        const size_t lda = m;
        MatrixND<T> result(m, lda);
        result.leftCols(k) = working;
        auto* a = reinterpret_cast<Tm*>(result.data());
        auto* tau = reinterpret_cast<const Tm*>(taus.data());

        if constexpr (isComplex) {
            if constexpr (T::Prec == Float32)
                check_lapack(LAPACKE_cungqr_64(Layout, m, n, k, a, lda, tau));
            else
                check_lapack(LAPACKE_zungqr_64(Layout, m, n, k, a, lda, tau));
        }
        else {
            if constexpr (T::Prec == Float32)
                check_lapack(LAPACKE_sorgqr_64(Layout, m, n, k, a, lda, tau));
            else
                check_lapack(LAPACKE_dorgqr_64(Layout, m, n, k, a, lda, tau));
        }
        return result;
    }
}
