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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Conjugate.h"

namespace Physica::Core {
    template<class ScalarType>
    class IterateSolver {
    public:
        using VectorType = Vector<ScalarType>;
        using RealType = typename ScalarType::RealType;
    private:
        VectorType solution;
    public:
        IterateSolver() = default;
        IterateSolver(size_t order);
        IterateSolver(const IterateSolver&) = default;
        IterateSolver(IterateSolver&&) noexcept = default;
        ~IterateSolver() = default;
        /* Operators */
        IterateSolver& operator=(IterateSolver solver) noexcept;
        /* Operations */
        template<class MatrixType>
        const VectorType& solve(const RValueMatrix<MatrixType>& A, VectorType b);
        void resize(size_t size) { solution.resize(size); }
        /* Getters */
        [[nodiscard]] size_t getOrder() const noexcept { return solution.getLength(); }
        [[nodiscard]] const VectorType& getSolution() const noexcept { return solution; }
        /* Helpers */
        void swap(IterateSolver& solver) noexcept;
    };

    template<class ScalarType>
    IterateSolver<ScalarType>::IterateSolver(size_t order) : solution(order) {}

    template<class ScalarType>
    IterateSolver<ScalarType>& IterateSolver<ScalarType>::operator=(IterateSolver solver) noexcept {
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
    template<class ScalarType>
    template<class MatrixType>
    const typename IterateSolver<ScalarType>::VectorType& IterateSolver<ScalarType>::solve(const RValueMatrix<MatrixType>& A, VectorType b) {
        solution = std::move(b);
        VectorType residual = -solution;
        VectorType p = solution;
        solution = ScalarType::Zero();
        VectorType& x = solution; //rename
        VectorType temp;

        RealType squaredNorm = residual.squaredNorm();
        while(squaredNorm > std::numeric_limits<ScalarType>::epsilon()) {
            temp = A * p;
            const ScalarType step = squaredNorm / (p.conjugate() * temp);
            x += step * p;
            residual += step * temp;
            const RealType next_squaredNorm = residual.squaredNorm();
            const RealType beta = next_squaredNorm / squaredNorm;
            squaredNorm = next_squaredNorm;
            p = beta * p - residual;
        }
        return solution;
    }

    template<class ScalarType>
    void IterateSolver<ScalarType>::swap(IterateSolver& solver) noexcept {
        solution.swap(solver.solution);
    }

    template<class ScalarType>
    inline void swap(IterateSolver<ScalarType>& solver1, IterateSolver<ScalarType>& solver2) noexcept {
        solver1.swap(solver2);
    }
}
