/*
 * Copyright 2026 Weibo He.
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

#include "Physica/Core/Exception/BadConvergenceException.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Givens.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/Orthogonalize.h"

namespace Physica {
    /**
     * Reference:
     * [1] SIAM J. Sci. Stat. Comput. 44, 856-869, (1986); https://doi.org/10.1137/0907058
     */
    template<Scalar T>
    class GMRES {
        using This = GMRES<T>;
        using Tr = T::RealType;
        constexpr static size_t Unlimited = std::numeric_limits<size_t>::max();
    private:
        MatrixND<T> krylov;
        MatrixND<T> hess;
        VectorND<T> residual;
        VectorND<T> coeffs;
        VectorND<T> buffer;
        size_t iteration{};
        size_t itelim = Unlimited;
        Tr res0;
        Tr res;
        Tr tolerance = std::numeric_limits<T>::epsilon();
    public:
        GMRES() = default;
        GMRES(size_t size, size_t krylovDim);
        GMRES(const This&) = default;
        GMRES(This&&) noexcept = default;
        ~GMRES() = default;
        /* Operators */
        This& operator=(This obj) noexcept;
        /* Operations */
        void solve(const Matrix auto& A, VectorND<T>& b);
        void solve(const Matrix auto& A, const VectorND<T>& b, VectorND<T>& x);

        void swap(This& obj) noexcept;
        /* Getters */
        [[nodiscard]] size_t getLength() const noexcept { return krylov.getRow(); }
        [[nodiscard]] size_t getKrylovDim() const noexcept { return krylov.getCol(); }
        [[nodiscard]] const auto& getResidual() noexcept { return residual; }
        [[nodiscard]] size_t getIteration() const noexcept { return iteration; }
        [[nodiscard]] size_t getMaxIteration() const noexcept;
        [[nodiscard]] bool isConverged() const noexcept;
        /* Setters */
        void setIterationLimit(size_t limit) noexcept { itelim = limit; }
        void setTolerance(Tr tol) noexcept { tolerance = tol; }
    private:
        void solve_impl(const Matrix auto& A, const VectorND<T>& b, VectorND<T>& x);
        [[nodiscard]] size_t spanKrylov(const Matrix auto& A);
    };

    template<Scalar T>
    GMRES<T>::GMRES(size_t size, size_t krylovDim)
            : krylov(size, krylovDim)
            , hess(krylovDim, krylovDim - 1, 0)
            , residual(size)
            , coeffs(krylovDim)
            , buffer(size) {
        assert(1 < krylovDim && krylovDim <= size);
    }

    template<Scalar T>
    auto GMRES<T>::operator=(This obj) noexcept -> This& {
        swap(obj);
        return *this;
    }

    template<Scalar T>
    void GMRES<T>::solve(const Matrix auto& A, VectorND<T>& b) {
        assert(b.getLength() == getLength());
        VectorND<T> x(getLength(), 0);
        b.assign(residual);
        res = res0 = b.norm();
        solve_impl(A, b, x);
        x.assign(b);
    }

    template<Scalar T>
    void GMRES<T>::solve(const Matrix auto& A, const VectorND<T>& b, VectorND<T>& x) {
        residual = b - A * x;
        res0 = b.norm();
        res = residual.norm();
        solve_impl(A, b, x);
    }

    template<Scalar T>
    void GMRES<T>::swap(This& obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        krylov.swap(obj.krylov);
        hess.swap(obj.hess);
        residual.swap(obj.residual);
        coeffs.swap(obj.coeffs);
        buffer.swap(obj.buffer);
        std::swap(iteration, obj.iteration);
        std::swap(itelim, obj.itelim);
        tolerance.swap(obj.tolerance);
    }

    template<Scalar T>
    size_t GMRES<T>::getMaxIteration() const noexcept {
        constexpr static int MaxIterationFactor = 8;
        return itelim == Unlimited ? std::max<size_t>(50, MaxIterationFactor * getLength()) : itelim;
    }

    template<Scalar T>
    bool GMRES<T>::isConverged() const noexcept {
        return res <= res0 * tolerance;
    }

    template<Scalar T>
    void GMRES<T>::solve_impl(const Matrix auto& A, const VectorND<T>& b, VectorND<T>& x) {
        assert(A.isSquare());
        assert(A.getRow() == krylov.getRow());
        const size_t maxIte = getMaxIteration();
        iteration = 0;
        while (!isConverged()) {
            if (iteration++ == maxIte) [[unlikely]]
                throw BadConvergenceException("[Error]: GMRES failed to converge");

            const size_t dimK = spanKrylov(A);
            const size_t dimV = dimK - (dimK == getKrylovDim());
            buffer.head(dimK).zeros();
            buffer[0] = res;
            for (size_t i = 0; i < dimV; ++i) {
                auto rotater = givens(hess.col(i), i, i + 1);
                applyGivens(rotater, hess.rightCols(i), i, i + 1);
                applyGivensCol(rotater, buffer, i, i + 1);
            }
            coeffs = hess.leftCols(dimV).triu().inv() * buffer.head(dimV);
            x += krylov.leftCols(dimV) * coeffs;
            residual = b - A * x;
            res = residual.norm();
        }
    }

    template<Scalar T>
    size_t GMRES<T>::spanKrylov(const Matrix auto& A) {
        krylov.col(0) = residual * reciprocal(res);
        for (size_t c = 0; c < hess.getCol(); ++c) {
            (A * krylov.col(c)).assign(buffer);

            auto v = krylov.col(c + 1);
            v = buffer;
            gram_schmidt(krylov.leftCols(c + 1), v, hess.col(c).head(c + 1));
            if (v.normInf() < buffer.normInf() * std::numeric_limits<T>::epsilon())
                return c + 1;
            v.toUnit();
            hess[c + 1, c] = v.conjugate() * buffer;
        }
        return getKrylovDim();
    }
}
