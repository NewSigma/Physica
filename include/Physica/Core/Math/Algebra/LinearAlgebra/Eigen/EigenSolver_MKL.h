/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "EigenSolver.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica {
    template<Scalar T, size_t Order>
    void EigenSolver<T, Order>::compute_mkl(const Matrix auto& source) {
        static_assert(T::Prec == Float32 || T::Prec == Float64);
        using Tm = decltype(std::declval<T>().toMKL());
        constexpr int Major = MatrixMajor::isColMatrix<decltype(source)>() ? MatrixMajor::Col : MatrixMajor::Row;
        constexpr int Layout = Major == MatrixMajor::Col ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
        using WorkingMatrixMKL = DenseMatrix<T, Major>;

        pre_compute(source);
        const size_t order = source.getRow();
        WorkingMatrixMKL working = source.values();
        auto* a = reinterpret_cast<Tm*>(working.data());
        if (getNeedEigenvectors()) {
            auto* vr = reinterpret_cast<Tm*>(rawEigenvectors.values().data());
            if constexpr (isComplex) {
                auto* w = reinterpret_cast<Tm*>(eigenvalues.data());
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgeev_64(Layout, 'N', 'V', order, a, order, w, nullptr, order, vr, order));
                else
                    check_lapack(LAPACKE_zgeev_64(Layout, 'N', 'V', order, a, order, w, nullptr, order, vr, order));
            }
            else {
                VectorND<T> ereal(order);
                VectorND<T> eimag(order);
                auto* wr = reinterpret_cast<Tm*>(ereal.data());
                auto* wi = reinterpret_cast<Tm*>(eimag.data());
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgeev_64(Layout, 'N', 'V', order, a, order, wr, wi, nullptr, order, vr, order));
                else
                    check_lapack(LAPACKE_dgeev_64(Layout, 'N', 'V', order, a, order, wr, wi, nullptr, order, vr, order));
                for (size_t i = 0; i < order; ++i)
                    eigenvalues[i] = Tc(ereal[i], eimag[i]);
            }

            if constexpr (Major == MatrixMajor::Row) {
                working  = rawEigenvectors.transpose();
                rawEigenvectors = working;
            }
        }
        else {
            MKL_INT64 sdim{};
            if constexpr (isComplex) {
                auto* w = reinterpret_cast<Tm*>(eigenvalues.data());
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_cgees_64(Layout, 'N', 'N', nullptr, order, a, order, &sdim, w, nullptr, order));
                else
                    check_lapack(LAPACKE_zgees_64(Layout, 'N', 'N', nullptr, order, a, order, &sdim, w, nullptr, order));
            }
            else {
                VectorND<T> ereal(order);
                VectorND<T> eimag(order);
                auto* wr = reinterpret_cast<Tm*>(ereal.data());
                auto* wi = reinterpret_cast<Tm*>(eimag.data());
                if constexpr (T::Prec == Float32)
                    check_lapack(LAPACKE_sgees_64(Layout, 'N', 'N', nullptr, order, a, order, &sdim, wr, wi, nullptr, order));
                else
                    check_lapack(LAPACKE_dgees_64(Layout, 'N', 'N', nullptr, order, a, order, &sdim, wr, wi, nullptr, order));
                for (size_t i = 0; i < order; ++i)
                    eigenvalues[i] = Tc(ereal[i], eimag[i]);
            }
        }
    }
}
