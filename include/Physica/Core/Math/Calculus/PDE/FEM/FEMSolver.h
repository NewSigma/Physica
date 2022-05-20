/*
 * Copyright 2022 WeiBo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/LinearEquations.h"

namespace Physica::Core {
    template<class ScalarType>
    class FEMSolver {
        using VectorType = Vector<ScalarType>;
        using MatrixType = DenseMatrix<ScalarType>;
    protected:
        MatrixType A;
        VectorType x;
        VectorType b;
    protected:
        FEMSolver() = default;
        FEMSolver(size_t n);
        /* Operations */
        void solve();
        void clear();
        void resize(size_t n);
    };

    template<class ScalarType>
    FEMSolver<ScalarType>::FEMSolver(size_t n) : A(n, n), x(n), b(n) {}

    template<class ScalarType>
    void FEMSolver<ScalarType>::solve() {
        const size_t n = x.getLength();
        MatrixType mat(n, n + 1);
        auto block = mat.leftCols(n);
        block = A;
        auto col = mat.col(n);
        col = b;
        LinearEquations equ(std::move(mat));
        equ.conjugateGradient();
        x = equ.getSolution();
    }

    template<class ScalarType>
    void FEMSolver<ScalarType>::clear() {
        A = ScalarType::Zero();
        x = ScalarType::Zero();
        b = ScalarType::Zero();
    }

    template<class ScalarType>
    void FEMSolver<ScalarType>::resize(size_t n) {
        A.resize(n, n);
        x.resize(n);
        b.resize(n);
    }
}
