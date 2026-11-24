/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../DenseLU.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica {
    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::compute_mkl() {
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr static int Layout = LAPACK_COL_MAJOR;
        const size_t m = getOrder();
        const size_t n = m;
        const size_t lda = m;
        auto* a = reinterpret_cast<Tm*>(working.data());
        int err{};
        if constexpr (Pivot) {
            Array<MKL_INT64> ipiv(m);
            if constexpr (T::isComplex()) {
                if constexpr (T::Prec == Float32)
                    err = LAPACKE_cgetrf_64(Layout, m, n, a, lda, ipiv.data());
                else
                    err = LAPACKE_zgetrf_64(Layout, m, n, a, lda, ipiv.data());
            }
            else {
                if constexpr (T::Prec == Float32)
                    err = LAPACKE_sgetrf_64(Layout, m, n, a, lda, ipiv.data());
                else
                    err = LAPACKE_dgetrf_64(Layout, m, n, a, lda, ipiv.data());
            }
            perm = PermType::fromLAPACK(ipiv);
        }
        else {
            if constexpr (T::isComplex()) {
                if constexpr (T::Prec == Float32)
                    err = LAPACKE_mkl_cgetrfnp_64(Layout, m, n, a, lda);
                else
                    err = LAPACKE_mkl_zgetrfnp_64(Layout, m, n, a, lda);
            }
            else {
                if constexpr (T::Prec == Float32)
                    err = LAPACKE_mkl_sgetrfnp_64(Layout, m, n, a, lda);
                else
                    err = LAPACKE_mkl_dgetrfnp_64(Layout, m, n, a, lda);
            }
        }

        // err > 0 implies a singular matrix, we do not care about it
        if (err < 0)
            check_lapack(err);
    }

    template<Scalar T, bool Pivot>
    void DenseLU<T, Pivot>::compute_mkl(const Matrix auto& source) {
        source.assign(working);
        compute_mkl();
    }
}
