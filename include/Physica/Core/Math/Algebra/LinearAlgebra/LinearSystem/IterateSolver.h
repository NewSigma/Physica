/*
 * Copyright 2022-2026 Weibo He.
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
        using Tr = T::RealType;
        constexpr static size_t Unlimited = std::numeric_limits<size_t>::max();
    private:
        VectorND<T> residual;
        VectorND<T> searchP;
        VectorND<T> dot;
        Tr squaredRes0;
        Tr squaredRes;
        Tr error = std::numeric_limits<T>::epsilon();
        size_t itelim = Unlimited;
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
        void cg(const Matrix auto& A, Vector auto&& b);
        void cg(std::invocable<const VectorND<T>&, VectorND<T>&> auto dotFunc, Vector auto&& b);
        void cgnr(const Matrix auto& A, Vector auto&& b);      
        void cgnr(std::invocable<const VectorND<T>&, VectorND<T>&> auto dotFunc, std::invocable<const VectorND<T>&, VectorND<T>&> auto dotTransFunc, Vector auto&& b);

        void resize(size_t size);
        void swap(This& __restrict obj) noexcept;
        /* Setters */
        void setError(Tr error_) noexcept { error = std::move(error_); }
        void setIterationLimit(size_t limit) noexcept { itelim = limit; }
    private:
        [[nodiscard]] bool isConverged() const noexcept;
        [[nodiscard]] size_t getMaxIterationCG() const noexcept;
    };

    template<Scalar T>
    IterateSolver<T>::IterateSolver(size_t size) {
        resize(size);
    }

    template<Scalar T>
    void IterateSolver<T>::cg(const Matrix auto& A, Vector auto&& b) {
        assert(A.isSquare());
        assert(A.getRow() == b.getLength());
        cg([&A](const VectorND<T>& v, VectorND<T>& dot) { dot = A * v; }, std::forward<decltype(b)>(b));
    }
    /**
     * Conjugate gradient(CG)
     * Pertinent for large-scale problem, \param dotFunc must be symmetric and positive definite
     * 
     * Reference:
     * [1] Nocedal J, Wright S J, Mikosch T V, et al. Numerical Optimization. Springer, 2006:112
     */
    template<Scalar T>
    void IterateSolver<T>::cg(std::invocable<const VectorND<T>&, VectorND<T>&> auto dotFunc, Vector auto&& b) {
        size_t length = b.getLength();
        if (dot.getLength() != length)
            resize(length);
        residual = -b;
        searchP = b;
        auto& x = b;
        x.zeros();

        squaredRes0 = squaredRes = residual.squaredNorm();
        const size_t maxite = getMaxIterationCG();
        size_t iteration = 0;
        while(!isConverged()) {
            dotFunc(searchP, dot);
            const Tr resA = (searchP.conjugate() * dot).real();
            assert((resA.isPositive() || !mustConverge) && "[Error]: The matrix is not positive definite");
            const Tr step = squaredRes / resA;
            x += step * searchP;
            residual += step * dot;

            const Tr next_squaredRes = residual.squaredNorm();
            const Tr beta = next_squaredRes / squaredRes;
            squaredRes = next_squaredRes;
            searchP = beta * searchP - residual;
            if (++iteration > maxite)
                break;
        }

        if (mustConverge && !isConverged()) [[unlikely]]
            throw BadConvergenceException("[Error]: CG cannot converge in the given iterations");
    }

    template<Scalar T>
    void IterateSolver<T>::cgnr(const Matrix auto& A, Vector auto&& b) {
        assert(A.isSquare());
        assert(A.getRow() == b.getLength());
        auto dotFunc = [&A](const VectorND<T>& v, VectorND<T>& dot) { (A * v).assign(dot); };
        auto dotTransFunc = [&A](const VectorND<T>& v, VectorND<T>& dot) { (A.hermite() * v).assign(dot); };
        cgnr(dotFunc, dotTransFunc, std::forward<decltype(b)>);
    }
    /**
     * Conjugate gradient normal equation residual(CGNR)
     * Applies to unsymmetric matrix, converges slow
     * 
     * References:
     * [1] Gene H. Golub, Charles F. Van Loan. Matrix computations 4th edition[M]. John Hopkins University Press, 2013:636-637
     */
    template<Scalar T>
    void IterateSolver<T>::cgnr(
            std::invocable<const VectorND<T>&, VectorND<T>&> auto dotFunc,
            std::invocable<const VectorND<T>&, VectorND<T>&> auto dotTransFunc,
            Vector auto&& b) {
        size_t length = b.getLength();
        if (dot.getLength() != length)
            resize(length);
        residual = -b;
        dotTransFunc(b, searchP);
        auto& x = b;
        x.zeros();

        squaredRes0 = squaredRes = searchP.squaredNorm();
        const size_t maxite = getMaxIterationCG();
        size_t iteration = 0;
        while(!isConverged()) {
            dotFunc(searchP, dot);
            const Tr step = squaredRes / dot.squaredNorm();
            x += step * searchP;
            residual += step * dot;

            dotTransFunc(residual, dot);
            const Tr next_squaredRes = dot.squaredNorm();
            const Tr beta = next_squaredRes / squaredRes;
            squaredRes = next_squaredRes;
            searchP = beta * searchP - residual;
            if (++iteration > maxite)
                break;
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
        std::swap(itelim, obj.itelim);
        std::swap(mustConverge, obj.mustConverge);
    }

    template<Scalar T>
    bool IterateSolver<T>::isConverged() const noexcept {
        return squaredRes < squaredRes0 * error;
    }

    template<Scalar T>
    size_t IterateSolver<T>::getMaxIterationCG() const noexcept {
        constexpr static int MaxIterationFactor = 2;
        return itelim == Unlimited ? (MaxIterationFactor * dot.getLength()) : itelim;
    }
}
