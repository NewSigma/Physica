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

namespace Physica::Core {
    template<class MatrixType, class VectorType>
    class IterateSolver {
        static_assert(std::is_convertible<MatrixType&, RValueMatrix<MatrixType>&>::value, "[Error]: Type 'MatrixType' must be a matrix");
        static_assert(std::is_convertible<VectorType&, LValueVector<VectorType>&>::value, "[Error]: Type 'VectorType' must be a lvector");
        using ScalarType  = typename MatrixType::ScalarType;
    public:
        MatrixType A;
        VectorType b;
    public:
        IterateSolver() = default;
        IterateSolver(MatrixType A_, VectorType b_);
        IterateSolver(const IterateSolver&) = default;
        IterateSolver(IterateSolver&&) noexcept = default;
        ~IterateSolver() = default;
        /* Operators */
        IterateSolver& operator=(IterateSolver solver) noexcept;
        /* Operations */
        void solve();
        /* Getters */
        [[nodiscard]] const VectorType& getSolution() { return b; }
        /* Helpers */
        void swap(IterateSolver& solver) noexcept;
    };

    template<class MatrixType, class VectorType>
    IterateSolver<MatrixType, VectorType>::IterateSolver(MatrixType A_, VectorType b_) : A(std::move(A_)), b(std::move(b_)) {}

    template<class MatrixType, class VectorType>
    IterateSolver<MatrixType, VectorType>& IterateSolver<MatrixType, VectorType>::operator=(IterateSolver solver) noexcept {
        swap(solver);
        return *this;
    }
    /**
     * Conjagate gradient method
     * Pertinent for large-scale problem, coefficient matrix must be symmetric and positive definite
     * 
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006.112
     */
    template<class MatrixType, class VectorType>
    void IterateSolver<MatrixType, VectorType>::solve() {
        VectorType residual = -b;
        VectorType p = b;
        b = ScalarType::Zero();
        VectorType& x = b; //rename
        VectorType temp;

        ScalarType squaredNorm = residual.squaredNorm();
        while(squaredNorm > std::numeric_limits<ScalarType>::epsilon()) {
            temp = A * p;
            const ScalarType step = squaredNorm / (p * temp);
            x += step * p;
            residual += step * temp;
            const ScalarType next_squaredNorm = residual.squaredNorm();
            const ScalarType beta = next_squaredNorm / squaredNorm;
            squaredNorm = next_squaredNorm;
            p = beta * p - residual;
        }
    }

    template<class MatrixType, class VectorType>
    void IterateSolver<MatrixType, VectorType>::swap(IterateSolver& solver) noexcept {
        A.swap(solver.A);
        b.swap(solver.b);
    }

    template<class MatrixType, class VectorType>
    inline void swap(IterateSolver<MatrixType, VectorType>& solver1,
                     IterateSolver<MatrixType, VectorType>& solver2) noexcept {
        solver1.swap(solver2);
    }
}
