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
    template<class ScalarType>
    class IterateSolver {
    public:
        using VectorType = Vector<ScalarType>; // IterateSolver is intended for large problem, so there is no need for fixed vectors
        using RealType = typename ScalarType::RealType;
    private:
        RealType error;
        size_t maxIteration;
    public:
        IterateSolver();
        IterateSolver(size_t maxIteration_, RealType error_);
        IterateSolver(const IterateSolver&) = default;
        IterateSolver(IterateSolver&&) noexcept = default;
        ~IterateSolver() = default;
        /* Operators */
        IterateSolver& operator=(IterateSolver solver) noexcept;
        /* Operations */
        template<class MatrixType>
        inline void solve(const RValueMatrix<MatrixType>& A, VectorType& b);
        template<class Functor, class AnyVector>
        void solve_functor(Functor dot_functor, LValueVector<AnyVector>& b);
        void resize([[maybe_unused]] size_t size) {}
        /* Helpers */
        void swap(IterateSolver& __restrict solver) noexcept;
        /* Setters */
        void setMaxIteration(size_t iteration) noexcept { maxIteration = iteration; }
    };

    template<class ScalarType>
    IterateSolver<ScalarType>::IterateSolver() : error(std::numeric_limits<ScalarType>::epsilon()), maxIteration(0) {}

    template<class ScalarType>
    IterateSolver<ScalarType>::IterateSolver(size_t maxIteration_, RealType error_)
            : error(std::move(error_)), maxIteration(maxIteration_) {}

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
    inline void IterateSolver<ScalarType>::solve(const RValueMatrix<MatrixType>& A, VectorType& b) {
        assert(A.getRow() == A.getColumn());
        assert(A.getRow() == b.getLength());
        solve_functor([&A](const VectorType& v, VectorType& dot) { dot = A.getDerived() * v; }, b);
    }

    template<class ScalarType>
    template<class Functor, class AnyVector>
    void IterateSolver<ScalarType>::solve_functor(Functor dot_functor, LValueVector<AnyVector>& b) {
        VectorType residual = -b.getDerived();
        VectorType p = b.getDerived();
        AnyVector& x = b.getDerived();
        x = ScalarType(0);
        VectorType dot(b.getLength());

        RealType squaredRes = residual.squaredNorm();
        const RealType resLimit = squaredRes * error;
        const bool noMaxIteration = maxIteration == 0;
        size_t iteration = 0;
        while((squaredRes > resLimit) && (iteration < maxIteration || noMaxIteration)) {
            dot_functor(p, dot);
            const RealType step = squaredRes / (p.conjugate() * dot).getReal();
            x += step * p;
            residual += step * dot;
            const RealType next_squaredRes = residual.squaredNorm();
            const RealType beta = next_squaredRes / squaredRes;
            squaredRes = next_squaredRes;
            p = beta * p - residual;

            ++iteration;
        }
    }

    template<class ScalarType>
    void IterateSolver<ScalarType>::swap(IterateSolver& __restrict solver) noexcept {
        assert(this != &solver && "[Error]: Self swap is likely a bug");
        error.swap(solver.error);
        std::swap(maxIteration, solver.maxIteration);
    }

    template<class ScalarType>
    inline void swap(IterateSolver<ScalarType>& __restrict solver1, IterateSolver<ScalarType>& __restrict solver2) noexcept {
        solver1.swap(solver2);
    }
}
