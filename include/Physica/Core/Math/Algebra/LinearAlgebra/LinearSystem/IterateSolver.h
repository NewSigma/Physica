/*
 * Copyright 2022-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Exception/BadConvergenceException.h"

namespace Physica {
    template<Scalar T>
    class IterateSolver {
        using This = IterateSolver<T>;
        using VectorType = VectorND<T>;
        using RealType = T::RealType;
    private:
        VectorType residual;
        VectorType searchP;
        VectorType dot;
        RealType squaredRes0;
        RealType squaredRes;
        RealType error = std::numeric_limits<T>::epsilon();
        size_t maxIteration = 0;
        size_t iteration;
    public:
        bool mustConverge = true;
    public:
        IterateSolver() = default;
        IterateSolver(size_t size);
        IterateSolver(const This&) = default;
        IterateSolver(This&&) noexcept = default;
        ~IterateSolver() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        inline void cg(const Matrix auto& A, VectorType& b);
        inline void cgnr(const Matrix auto& A, VectorType& b);

        void cg_functor(std::invocable<const VectorType&, VectorType&> auto dotFunc, Vector auto& b);
        void cgnr_functor(std::invocable<const VectorType&, VectorType&> auto dotFunc, std::invocable<const VectorType&, VectorType&> auto dotTransFunc, Vector auto& b);

        void resize(size_t size);
        void swap(This& __restrict obj) noexcept;
        /* Setters */
        void setError(RealType error_) noexcept { error = std::move(error_); }
        void setMaxIteration(size_t iteration) noexcept { maxIteration = iteration; }
    private:
        [[nodiscard]] bool isConverged() const { return squaredRes < squaredRes0 * error; }
        [[nodiscard]] bool shouldStop() const { return (maxIteration != 0) && (iteration > maxIteration); }
    };

    template<Scalar T>
    IterateSolver<T>::IterateSolver(size_t size) {
        resize(size);
    }

    template<Scalar T>
    inline void IterateSolver<T>::cg(const Matrix auto& A, VectorType& b) {
        assert(A.getRow() == A.getCol());
        assert(A.getRow() == b.getLength());
        cg_functor([&A](const VectorType& v, VectorType& dot) { dot = A * v; }, b);
    }

    template<Scalar T>
    inline void IterateSolver<T>::cgnr(const Matrix auto& A, VectorType& b) {
        assert(A.getRow() == A.getCol());
        assert(A.getRow() == b.getLength());
        auto dotFunc = [&A](const VectorType& v, VectorType& dot) { dot = A * v; };
        auto dotTransFunc = [&A](const VectorType& v, VectorType& dot) { dot = A.transpose() * v; };
        cgnr_functor(dotFunc, dotTransFunc, b);
    }
    /**
     * Conjagate gradient(CG)
     * Pertinent for large-scale problem, \param dotFunc must be symmetric and positive definite
     * 
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006:112
     */
    template<Scalar T>
    void IterateSolver<T>::cg_functor(std::invocable<const VectorType&, VectorType&> auto dotFunc, Vector auto& b) {
        if (dot.getLength() != b.getLength())
            resize(b.getLength());
        residual = -b;
        searchP = b;
        auto& x = b;
        x.zeros();

        squaredRes0 = squaredRes = residual.squaredNorm();
        iteration = 0;
        while(!isConverged() && !shouldStop()) {
            dotFunc(searchP, dot);
            const RealType step = squaredRes / (searchP.conjugate() * dot).real();
            x += step * searchP;
            residual += step * dot;

            const RealType next_squaredRes = residual.squaredNorm();
            const RealType beta = next_squaredRes / squaredRes;
            squaredRes = next_squaredRes;
            searchP = beta * searchP - residual;
            ++iteration;
        }

        if (mustConverge && !isConverged()) [[unlikely]]
            throw BadConvergenceException("[Error]: CG cannot converge in the given iterations");
    }
    /**
     * Conjagate gradient normal equation residual(CGNR)
     * Applies to unsymmetric matrix, converges slow
     * 
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:636-637
     */
    template<Scalar T>
    void IterateSolver<T>::cgnr_functor(
            std::invocable<const VectorType&, VectorType&> auto dotFunc,
            std::invocable<const VectorType&, VectorType&> auto dotTransFunc,
            Vector auto& b) {
        if (dot.getLength() != b.getLength())
            resize(b.getLength());
        residual = -b;
        dotTransFunc(b, searchP);
        auto& x = b;
        x.zeros();

        squaredRes0 = squaredRes = searchP.squaredNorm();
        iteration = 0;
        while(!isConverged() && !shouldStop()) {
            dotFunc(searchP, dot);
            const RealType step = squaredRes / dot.squaredNorm();
            x += step * searchP;
            residual += step * dot;

            dotTransFunc(residual, dot);
            const RealType next_squaredRes = dot.squaredNorm();
            const RealType beta = next_squaredRes / squaredRes;
            squaredRes = next_squaredRes;
            searchP = beta * searchP - residual;
            ++iteration;
        }

        if (mustConverge && !isConverged()) [[unlikely]]
            throw BadConvergenceException("[Error]: CGNR cannot converge in the given iterations");
    }

    template<Scalar T>
    void IterateSolver<T>::resize(size_t size) {
        residual.resize(size);
        searchP.resize(size);
        dot.resize(size);
    }

    template<Scalar T>
    void IterateSolver<T>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        residual.swap(obj.residual);
        searchP.swap(obj.searchP);
        dot.swap(obj.dot);
        squaredRes0.swap(obj.squaredRes0);
        squaredRes.swap(obj.squaredRes);
        error.swap(obj.error);
        std::swap(maxIteration, obj.maxIteration);
        std::swap(iteration, obj.iteration);
        std::swap(mustConverge, obj.mustConverge);
    }
}

namespace std {
    template<Physica::Scalar T>
    inline void swap(Physica::IterateSolver<T>& __restrict solver1,
                     Physica::IterateSolver<T>& __restrict solver2) noexcept {
        solver1.swap(solver2);
    }
}
