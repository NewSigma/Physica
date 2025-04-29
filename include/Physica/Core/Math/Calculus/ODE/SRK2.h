/*
 * Copyright 2022 Weibo He.
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
     * Defination of RandomFunc:
     * VectorType RandomFunc(T x, VectorType y);
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

        template<class Function, class RandomFunc>
        void solve(Function func, RandomFunc random);

        template<class Function, class RandomFunc>
        static inline void step(T stepSize, T& x, VectorType& sol, Function func, RandomFunc random);
    };

    template<Scalar T, size_t Dim>
    template<class Function, class RandomFunc>
    void SRK2<T, Dim>::solve(Function func, RandomFunc random) {
        const size_t col_1 = Base::solution.getCol() - 1;
        for (size_t i = 0; i < col_1; ++i) {
            T temp = Base::x[i];
            Base::solution.asArray()[i + 1] = Base::solution.col(i);
            step(Base::stepSize, temp, Base::solution.asArray()[i + 1], func, random);
            Base::x[i + 1] = temp;
        }
    }

    template<Scalar T, size_t Dim>
    template<class Function, class RandomFunc>
    inline void SRK2<T, Dim>::step(T stepSize, T& x, VectorType& sol, Function func, RandomFunc random) {
        using FunctionResult = std::invoke_result<Function, T, VectorType>::type;
        using RandFunctionResult = std::invoke_result<RandomFunc, T, VectorType>::type;
        static_assert(FunctionResult::SizeAtCompile == Dim, "[Possible optimization]: Dimention between ODESolver and functor do not match");
        static_assert(RandFunctionResult::SizeAtCompile == Dim, "[Possible optimization]: Dimention between ODESolver and functor do not match");
        const VectorType randVec = random(x, sol);
        VectorType term1 = func(x, sol);
        x += stepSize;
        VectorType term2 = func(x, sol + stepSize * term1 + randVec);
        sol += (term1 + term2) * (stepSize * T(0.5)) + randVec;
    }
}
