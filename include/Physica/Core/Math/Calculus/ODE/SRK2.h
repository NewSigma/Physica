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

#include "ODESolver.h"

namespace Physica::Core {
    /**
     * Apply to wight noise only
     * Defination of RandomFunc:
     * VectorType RandomFunc(T x, VectorType y);
     * 
     * Reference:
     * [1] R. L. Honeycutt, Stochastic Runge-Kutta algorithm: I. White noise, Phys. Rev. A 45, 600 (1992).
     */
    template<class ScalarType>
    class SRK2 : public ODESolver<ScalarType> {
        using Base = ODESolver<ScalarType>;
        using typename Base::VectorType;
    public:
        using Base::Base;
        ~SRK2() = default;

        template<class Function, class RandomFunc>
        void solve(Function func, RandomFunc random);

        template<class Function, class RandomFunc>
        static void step(ScalarType stepSize, ScalarType& x, VectorType& sol, Function func, RandomFunc random);
    };

    template<class ScalarType>
    template<class Function, class RandomFunc>
    void SRK2<ScalarType>::solve(Function func, RandomFunc random) {
        assert(Base::x.getLength() == 1);
        const size_t column_1 = Base::solution.getColumn() - 1;
        for (size_t i = 0; i < column_1; ++i) {
            ScalarType temp = Base::x[i];
            Base::solution[i + 1] = Base::solution[i];
            step(Base::stepSize, temp, Base::solution[i + 1], func, random);
            Base::x.append(temp);
        }
    }

    template<class ScalarType>
    template<class Function, class RandomFunc>
    void SRK2<ScalarType>::step(ScalarType stepSize, ScalarType& x, VectorType& sol, Function func, RandomFunc random) {
        const VectorType randVec = random(x, sol);
        VectorType term1 = func(x, sol);
        x += stepSize;
        VectorType term2 = func(x, sol + stepSize * term1 + randVec);
        sol += (term1 + term2) * (stepSize * ScalarType(0.5)) + randVec;
    }
}
