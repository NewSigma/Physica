/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Inverse.h"

namespace Physica {
    template<Scalar T, bool Pivot>
    void Inverse<DenseLU<T, Pivot>>::assign_mkl(Matrix auto& target) const {
        static_assert(target.getMinorStrideAtCompile() == 1, "[Error]: LAPACK expects a contiguous minor direction");
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr int Layout = LAPACK_COL_MAJOR;
        target.assert_assign(*this);
        lu.getMatrixLU().assign(target);

        const size_t n = getOrder();
        auto* a = reinterpret_cast<Tm*>(target.data_handle());
        const size_t lda = target.getMajorStride();
        if constexpr (Pivot) {
            auto buffer = lu.getPerm().template toLAPACK<MKL_INT64>();
            const auto* ipiv = buffer.data();
            if constexpr (isComplex()) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgetri_64(Layout, n, a, lda, ipiv));
                else
                    check_lapack(LAPACKE_zgetri_64(Layout, n, a, lda, ipiv));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgetri_64(Layout, n, a, lda, ipiv));
                else
                    check_lapack(LAPACKE_dgetri_64(Layout, n, a, lda, ipiv));
            }
        }
        else {
            if constexpr (isComplex()) {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_mkl_cgetrinp_64(Layout, n, a, lda));
                else
                    check_lapack(LAPACKE_mkl_zgetrinp_64(Layout, n, a, lda));
            }
            else {
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_mkl_sgetrinp_64(Layout, n, a, lda));
                else
                    check_lapack(LAPACKE_mkl_dgetrinp_64(Layout, n, a, lda));
            }
        }
    }
}
