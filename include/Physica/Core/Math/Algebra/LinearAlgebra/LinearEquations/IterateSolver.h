/*
 * Copyright 2022-2024 Weibo He.
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

#include <iostream>

namespace Physica::Core {
    template<class ScalarType>
    class IterateSolver {
        using This = IterateSolver<ScalarType>;
        using VectorType = Vector<ScalarType>;
        using RealType = typename ScalarType::RealType;
    private:
        VectorType residual;
        VectorType searchP;
        VectorType dot;
        RealType squaredRes0;
        RealType squaredRes;
        RealType error;
        size_t maxIteration;
        size_t iteration;
    public:
        IterateSolver();
        IterateSolver(size_t maxIteration_, RealType error_);
        IterateSolver(size_t size, size_t maxIteration_, RealType error_);
        IterateSolver(const This&) = default;
        IterateSolver(This&&) noexcept = default;
        ~IterateSolver() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        /* Operations */
        template<class MatrixType>
        inline void cg(const RValueMatrix<MatrixType>& A, VectorType& b);
        template<class MatrixType>
        inline void cgnr(const RValueMatrix<MatrixType>& A, VectorType& b);

        template<class Functor, class LVector>
        void cg_functor(Functor dotFunc, LVector& b);
        template<class Functor1, class Functor2, class LVector>
        void cgnr_functor(Functor1 dotFunc, Functor2 dotTransFunc, LVector& b);

        void resize(size_t size);
        void swap(This& __restrict obj) noexcept;
        /* Setters */
        void setMaxIteration(size_t iteration) noexcept { maxIteration = iteration; }
    private:
        [[nodiscard]] bool isConverged() const { return squaredRes < squaredRes0 * error; }
        [[nodiscard]] bool shouldStop() const { return (maxIteration != 0) && (iteration > maxIteration); }
    };

    template<class ScalarType>
    IterateSolver<ScalarType>::IterateSolver() : IterateSolver(0, std::numeric_limits<ScalarType>::epsilon()) {}

    template<class ScalarType>
    IterateSolver<ScalarType>::IterateSolver(size_t maxIteration_, RealType error_)
            : error(std::move(error_))
            , maxIteration(maxIteration_) {}

    template<class ScalarType>
    IterateSolver<ScalarType>::IterateSolver(size_t size, size_t maxIteration_, RealType error_)
            : IterateSolver(maxIteration_, error_) {
        resize(size);
    }

    template<class ScalarType>
    template<class MatrixType>
    inline void IterateSolver<ScalarType>::cg(const RValueMatrix<MatrixType>& A, VectorType& b) {
        assert(A.getRow() == A.getColumn());
        assert(A.getRow() == b.getLength());
        cg_functor([&A](const VectorType& v, VectorType& dot) { dot = A.getDerived() * v; }, b);
    }

    template<class ScalarType>
    template<class MatrixType>
    inline void IterateSolver<ScalarType>::cgnr(const RValueMatrix<MatrixType>& A, VectorType& b) {
        assert(A.getRow() == A.getColumn());
        assert(A.getRow() == b.getLength());
        auto dotFunc = [&A](const VectorType& v, VectorType& dot) { dot = A.getDerived() * v; };
        auto dotTransFunc = [&A](const VectorType& v, VectorType& dot) { dot = A.getDerived().transpose() * v; };
        cgnr_functor(dotFunc, dotTransFunc, b);
    }
    /**
     * Conjagate gradient(CG)
     * Pertinent for large-scale problem, \param dotFunc must be symmetric and positive definite
     * 
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006.112
     */
    template<class ScalarType>
    template<class Functor, class LVector>
    void IterateSolver<ScalarType>::cg_functor(Functor dotFunc, LVector& b) {
        static_assert(std::is_base_of<LValueVector<LVector>, LVector>::value, "[Error]: This is not a LValueVector");
        residual = -b;
        searchP = b;
        LVector& x = b;
        x = ScalarType(0);

        squaredRes0 = squaredRes = residual.squaredNorm();
        iteration = 0;
        while(!isConverged() && !shouldStop()) {
            dotFunc(searchP, dot);
            const RealType step = squaredRes / (searchP.conjugate() * dot).getReal();
            x += step * searchP;
            residual += step * dot;

            const RealType next_squaredRes = residual.squaredNorm();
            const RealType beta = next_squaredRes / squaredRes;
            squaredRes = next_squaredRes;
            searchP = beta * searchP - residual;
            ++iteration;
        }

        if (!isConverged()) [[unlikely]]
            throw BadConvergenceException("[Error]: CG cannot converge in the given iterations");
    }
    /**
     * Conjagate gradient normal equation residual(CGNR)
     * Applies to unsymmetric matrix, converges slow
     * 
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:636-637
     */
    template<class ScalarType>
    template<class Functor1, class Functor2, class LVector>
    void IterateSolver<ScalarType>::cgnr_functor(Functor1 dotFunc, Functor2 dotTransFunc, LVector& b) {
        static_assert(std::is_base_of<LValueVector<LVector>, LVector>::value, "[Error]: This is not a LValueVector");
        residual = -b;
        dotTransFunc(b, searchP);
        LVector& x = b;
        x = ScalarType(0);

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

        if (!isConverged()) [[unlikely]]
            throw BadConvergenceException("[Error]: CGNR cannot converge in the given iterations");
    }

    template<class ScalarType>
    void IterateSolver<ScalarType>::resize(size_t size) {
        residual.resize(size);
        searchP.resize(size);
        dot.resize(size);
    }

    template<class ScalarType>
    void IterateSolver<ScalarType>::swap(This& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        residual.swap(obj.residual);
        searchP.swap(obj.searchP);
        dot.swap(obj.dot);
        squaredRes0.swap(obj.squaredRes0);
        squaredRes.swap(obj.squaredRes);
        error.swap(obj.error);
        std::swap(maxIteration, obj.maxIteration);
        std::swap(iteration, obj.iteration);
    }
}

namespace std {
    template<class ScalarType>
    inline void swap(Physica::Core::IterateSolver<ScalarType>& __restrict solver1,
                     Physica::Core::IterateSolver<ScalarType>& __restrict solver2) noexcept {
        solver1.swap(solver2);
    }
}
