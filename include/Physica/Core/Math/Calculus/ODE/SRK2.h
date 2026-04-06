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

#include "ODESolver.h"

namespace Physica {
    /**
     * Apply to wight noise only
     * 
     * Reference:
     * [1] Phys. Rev. A 45, 600 (1992); https://doi.org/10.1103/PhysRevA.45.600
     */
    template<Scalar T, size_t Dim>
    class SRK2 : public ODESolver<T, Dim> {
        using Base = ODESolver<T, Dim>;
    public:
        using typename Base::VectorType;
    public:
        using Base::Base;
        ~SRK2() = default;
        /* Operations */
        void solve(std::invocable<T, VectorType> auto fn, std::invocable<T, VectorType> auto  random);
        /* Static members */
        static void step(T stepSize, T& x, Vector auto&& sol, std::invocable<T, VectorType> auto fn, std::invocable<T, VectorType> auto random);
    };

    template<Scalar T, size_t Dim>
    void SRK2<T, Dim>::solve(std::invocable<T, VectorType> auto fn, std::invocable<T, VectorType> auto random) {
        const size_t col_1 = Base::solution.getCol() - 1;
        for (size_t i = 0; i < col_1; ++i) {
            T temp = Base::x[i];
            Base::solution.col(i + 1) = Base::solution.col(i);
            step(Base::stepSize, temp, Base::solution.col(i + 1), fn, random);
            Base::x[i + 1] = temp;
        }
    }

    template<Scalar T, size_t Dim>
    void SRK2<T, Dim>::step(
            T stepSize,
            T& x,
            Vector auto&& sol,
            std::invocable<T, VectorType> auto fn,
            std::invocable<T, VectorType> auto random) {
        using FunctionResult = std::invoke_result<decltype(fn), T, VectorType>::type;
        using RandFunctionResult = std::invoke_result<decltype(random), T, VectorType>::type;
        static_assert(FunctionResult::SizeAtCompile == Dim, "[Possible optimization]: Dimention between ODESolver and functor do not match");
        static_assert(RandFunctionResult::SizeAtCompile == Dim, "[Possible optimization]: Dimention between ODESolver and functor do not match");
        const VectorType randVec = random(x, sol);
        VectorType term1 = fn(x, sol);
        x += stepSize;
        VectorType term2 = fn(x, sol + stepSize * term1 + randVec);
        sol += (term1 + term2) * (stepSize * T(0.5)) + randVec;
    }
}
