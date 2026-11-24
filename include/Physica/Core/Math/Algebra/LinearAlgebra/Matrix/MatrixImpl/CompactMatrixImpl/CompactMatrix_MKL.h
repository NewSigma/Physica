/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <mkl_lapacke.h>
#include <mkl_trans.h>
#include "../CompactMatrix.h"

namespace Physica {
    template<class Derived>
    void CompactMatrix<Derived>::assign_mkl(Matrix auto&& target) const noexcept {
        static_assert(std::remove_cvref_t<decltype(target)>::isCompact(), "[Error]: MKL expects compact matrix");
        static_assert(!MatrixMajor::isSameMajor<This, decltype(target)>(), "[Error]: Do not need transpose, use assign_base() instead");
        using Tm = decltype(std::declval<T>().toMKL());
        target.assert_assign_lapack(Base::getDerived());

        constexpr char ordering = 'R'; // Ordering does not matter
        constexpr char trans = 'T';
        const size_t rows = Base::getRow();
        const size_t cols = Base::getCol();
        const auto* A = reinterpret_cast<const Tm*>(data_handle());
        auto* B = reinterpret_cast<Tm*>(target.data_handle());
        const size_t lda = Base::getMajorStride();
        const size_t ldb = target.getMajorStride();
        if constexpr (T::isComplex()) {
            const auto alpha = T(1).toMKL();
            if constexpr (T::Prec == Float32)
                mkl_comatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
            else
                mkl_zomatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
        }
        else {
            constexpr double alpha = 1;
            if constexpr (T::Prec == Float32)
                mkl_somatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
            else
                mkl_domatcopy(ordering, trans, rows, cols, alpha, A, lda, B, ldb);
        }
    }

    template<class Derived>
    auto CompactMatrix<Derived>::balance_mkl() -> VectorND<T> {
        using Tm = decltype(std::declval<T>().toMKL());
        Base::assert_balance();

        constexpr int Layout = Derived::getMajor() == MatrixMajor::Row ? CblasRowMajor : CblasColMajor;
        constexpr char job = 'S';
        size_t n = Base::getRow();
        auto* a = reinterpret_cast<Tm*>(data_handle());
        size_t lda = Base::getMajorStride();
        lapack_int ilo{}, ihi{};
        VectorND<T> result(n);
        auto* scale = reinterpret_cast<Tm*>(result.data());
        if constexpr (isComplex()) {
            if constexpr (T::Prec == Float32)
                LAPACKE_cgebal(Layout, job, n, a, lda, &ilo, &ihi, scale);
            else
                LAPACKE_zgebal(Layout, job, n, a, lda, &ilo, &ihi, scale);
        }
        else {
            if constexpr (T::Prec == Float32)
                LAPACKE_sgebal(Layout, job, n, a, lda, &ilo, &ihi, scale);
            else
                LAPACKE_dgebal(Layout, job, n, a, lda, &ilo, &ihi, scale);
        }
        return reciprocal(result);
    }
}
