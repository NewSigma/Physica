/*
 * Copyright 2025 Weibo He.
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

#include "../SparseLU.h"
#include "Physica/Core/Exception/MKL/Pardiso.h"

namespace Physica {
    template<Scalar T>
    void SparseLU<T>::compute(const Matrix auto& source) {
        spmat = SparseMatrix<T>(source);
        runPardiso(Phase::AnalyseFactor);
    }

    template<Scalar T>
    void SparseLU<T>::analyse(const Matrix auto& source) {
        spmat = SparseMatrix<T>(source);
        runPardiso(Phase::Analyse);
    }

    template<Scalar T>
    void SparseLU<T>::factorize(const Matrix auto& source) {
        spmat = SparseMatrix<T>(source);
        runPardiso(Phase::Factorize);
    }

    template<Scalar T>
    MatrixND<T> SparseLU<T>::solve_mkl(const MatrixND<T>& rhs) {
        assert(rhs.getRow() == getOrder());
        MatrixND<T> result(rhs.getRow(), rhs.getCol());
        runPardiso(Phase::Solve, rhs.getCol(), const_cast<T*>(rhs.data()), result.data());
        return result;
    }

    template<Scalar T>
    void SparseLU<T>::initialize(bool needDiag) noexcept {
        void* pt = static_cast<void*>(handles.data());
        const int64_t mtype = getMatrixType();
        auto* iparm = (MKL_INT64*)params.data();
        pardisoinit(pt, (MKL_INT64*)&mtype, iparm);
        params[34] = 1;
        if (needDiag) {
            params[10] = 0;
            params[55] = 1;
        }
    }

    template<Scalar T>
    void SparseLU<T>::runPardiso(Phase phase, int64_t nrhs, T* b, T* x) {
        assert(((nrhs == 0) == (b == nullptr)) && "[Error]: Params inconsistent");
        assert(!isSolvePhase(phase) || nrhs > 0 && "[Error]: Empty RHS");
        void* pt = static_cast<void*>(handles.data());
        const MKL_INT64 maxfct = 1;
        const MKL_INT64 mtype = getMatrixType();
        const MKL_INT64 n = getOrder();
        const void* a = spmat.getElements().data();
        const auto* ia = (MKL_INT64*)spmat.getMajorStarts().data();
        const auto* ja = (MKL_INT64*)spmat.getMinorIndexes().data();
        MKL_INT64* iperm = nullptr;
        auto* iparm = (MKL_INT64*)params.data();
        const MKL_INT64 msglvl = 0; // We do not care about logs
        MKL_INT64 err{};
        pardiso_64(pt, &maxfct, &maxfct, &mtype, (MKL_INT64*)&phase, &n, a, ia, ja, iperm, reinterpret_cast<MKL_INT64*>(&nrhs), iparm, &msglvl, b, x, &err);
        check_pardiso(static_cast<int>(err));

        if (getNeedDiag() && isFactorizePhase(phase)) {
            pardiso_getdiag(pt, diags.data(), buffer.data(), &maxfct, &err);
            assert(err == 0 && "[Error]: We should have set iparm[55]");
        }
    }

    template<Scalar T>
    constexpr bool SparseLU<T>::isFactorizePhase(Phase p) noexcept {
        switch (p) {
        case Analyse:
            return false;
        case Factorize:
        case AnalyseFactorize:
        case AnalyseFactorizeSolve:
            return true;
        case Solve:
        case Destroy:
            return false;
        default:
            unreachable();
        }
    }

    template<Scalar T>
    constexpr bool SparseLU<T>::isSolvePhase(Phase p) noexcept {
        switch (p) {
        case Analyse:
        case Factorize:
        case AnalyseFactorize:
            return false;
        case AnalyseFactorizeSolve:
        case Solve:
            return true;
        case Destroy:
            return false;
        default:
            unreachable();
        }
    }

    template<Scalar T>
    constexpr int SparseLU<T>::getMatrixType() noexcept {
        return T::isComplex ? 13 : 11;
    }
}
