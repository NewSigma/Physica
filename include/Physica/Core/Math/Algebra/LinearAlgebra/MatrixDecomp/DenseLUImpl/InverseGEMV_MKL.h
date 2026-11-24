/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "InverseGEMV.h"

namespace Physica {
    template<Matrix M, Vector V> requires(instanceof_tx<M, Inverse> && instanceof_tx<typename Traits<M>::ExprType, DenseLU>)
    void GEMV<M, V>::assign_mkl(Vector auto& target) const {
        using Tm = decltype(std::declval<T>().toMKL());
        v.assign(target);

        constexpr int Layout = LAPACK_COL_MAJOR;
        constexpr char Trans = 'N';
        const auto& lu = m.getDenseLU();
        size_t n = getLength();
        const auto* a = reinterpret_cast<const Tm*>(lu.getMatrixLU().data());
        auto* b = reinterpret_cast<Tm*>(target.data());
        if constexpr (lu.isPivot()) {
            auto buffer = lu.getPerm().template toLAPACK<MKL_INT64>();
            const auto* ipiv = buffer.data();
            if constexpr (Base::isComplex()) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgetrs_work_64(Layout, Trans, n, 1, a, n, ipiv, b, n));
                else
                    check_lapack(LAPACKE_zgetrs_work_64(Layout, Trans, n, 1, a, n, ipiv, b, n));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgetrs_work_64(Layout, Trans, n, 1, a, n, ipiv, b, n));
                else
                    check_lapack(LAPACKE_dgetrs_work_64(Layout, Trans, n, 1, a, n, ipiv, b, n));
            }
        }
        else {
            Array<MKL_INT64> perm(n);
            for (size_t i = 0; i < n; ++i)
                perm[i] = i + 1;
            auto* ipiv = perm.data();
            if constexpr (Base::isComplex()) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgetrs_work_64(Layout, Trans, n, 1, a, n, ipiv, b, n));
                else
                    check_lapack(LAPACKE_zgetrs_work_64(Layout, Trans, n, 1, a, n, ipiv, b, n));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgetrs_work_64(Layout, Trans, n, 1, a, n, ipiv, b, n));
                else
                    check_lapack(LAPACKE_dgetrs_work_64(Layout, Trans, n, 1, a, n, ipiv, b, n));
            }
        }
    }
}
