/*
 * Copyright 2024 Weibo He.
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

#include "EigenSolver.h"
#include "Physica/Core/Exception/MKL/Lapack.h"

namespace Physica::Core {
    template<Scalar T, size_t Order>
    template<Matrix M>
    void EigenSolver<T, Order>::compute_mkl(const M& source, bool computeEigenvectors_) {
        static_assert(T::Option == Float32 || T::Option == Float64);
        constexpr int Major = MatrixOption::isRowMatrix<M>() ? MatrixOption::Row : MatrixOption::Col;
        constexpr int Layout = Major == MatrixOption::Row ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
        using WorkingMatrixMKL = DenseMatrix<T, Major | MatrixOption::Element>;
        using Tm = std::conditional<isComplex, typename ComplexType::MKL_Complex, typename T::MachineType>::type;

        pre_compute(source, computeEigenvectors_);
        const size_t order = source.getRow();
        WorkingMatrixMKL working = source;
        auto* a = reinterpret_cast<Tm*>(working.data());
        if (computeEigenvectors) {
            auto* vl = reinterpret_cast<Tm*>(rawEigenvectors.data());
            if constexpr (isComplex) {
                auto* w = reinterpret_cast<Tm*>(eigenvalues.data());
                if constexpr (T::Option == Float32)
                    check_lapack(LAPACKE_cgeev_64(Layout, 'N', 'V', order, a, order, w, nullptr, order, vl, order));
                else
                    check_lapack(LAPACKE_zgeev_64(Layout, 'N', 'V', order, a, order, w, nullptr, order, vl, order));
            }
            else {
                VectorND<T> ereal(order);
                VectorND<T> eimag(order);
                auto* wr = reinterpret_cast<Tm*>(ereal.data());
                auto* wi = reinterpret_cast<Tm*>(eimag.data());
                if constexpr (T::Option == Float32)
                    check_lapack(LAPACKE_sgeev_64(Layout, 'N', 'V', order, a, order, wr, wi, nullptr, order, vl, order));
                else
                    check_lapack(LAPACKE_dgeev_64(Layout, 'N', 'V', order, a, order, wr, wi, nullptr, order, vl, order));
                for (size_t i = 0; i < order; ++i)
                    eigenvalues[i] = ComplexType(ereal[i], eimag[i]);
            }
        }
        else {
            MKL_INT64 sdim;
            if constexpr (isComplex) {
                auto* w = reinterpret_cast<Tm*>(eigenvalues.data());
                if constexpr (T::Option == Float32)
                    check_lapack(LAPACKE_cgees_64(Layout, 'N', 'N', nullptr, order, a, order, &sdim, w, nullptr, order));
                else
                    check_lapack(LAPACKE_zgees_64(Layout, 'N', 'N', nullptr, order, a, order, &sdim, w, nullptr, order));
            }
            else {
                VectorND<T> ereal(order);
                VectorND<T> eimag(order);
                auto* wr = reinterpret_cast<Tm*>(ereal.data());
                auto* wi = reinterpret_cast<Tm*>(eimag.data());
                if constexpr (T::Option == Float32)
                    check_lapack(LAPACKE_sgees_64(Layout, 'N', 'N', nullptr, order, a, order, &sdim, wr, wi, nullptr, order));
                else
                    check_lapack(LAPACKE_dgees_64(Layout, 'N', 'N', nullptr, order, a, order, &sdim, wr, wi, nullptr, order));
                for (size_t i = 0; i < order; ++i)
                    eigenvalues[i] = ComplexType(ereal[i], eimag[i]);
            }
        }
    }
}
