/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_INTEGRATESOLVER_H
#define PHYSICA_INTEGRATESOLVER_H

#include "Integrate.h"

namespace Physica::Core {
    enum IntegrateMethod {
        Rectangular,
        Ladder,
        Simpson,
        Tanh_Sinh
    };

    template<IntegrateMethod method, size_t dim, ScalarOption option = MultiPrecision, bool errorTrack = true>
    class IntegrateSolver;
    //////////////////////////////////Rectangular//////////////////////////////////
    template<size_t dim, ScalarOption option, bool errorTrack>
    class IntegrateSolver<Rectangular, dim, option, errorTrack> {

    };

    template<ScalarOption option, bool errorTrack>
    class IntegrateSolver<Rectangular, 1, option, errorTrack> {
        Scalar<option, false> stepSize;
    public:
        explicit IntegrateSolver(Scalar<option, false> stepSize);
        /* Operations */
        Scalar<option, errorTrack> solve(const Integrate<1, option, errorTrack>& i) const;
    };
    //////////////////////////////////Ladder//////////////////////////////////
    template<size_t dim, ScalarOption option, bool errorTrack>
    class IntegrateSolver<Ladder, dim, option, errorTrack> {

    };

    template<ScalarOption option, bool errorTrack>
    class IntegrateSolver<Ladder, 1, option, errorTrack> {
        Scalar<option, false> stepSize;
    public:
        explicit IntegrateSolver(Scalar<option, false> stepSize);
        /* Operations */
        Scalar<option, errorTrack> solve(const Integrate<1, option, errorTrack>& i) const;
    };
    //////////////////////////////////Simpson//////////////////////////////////
    template<size_t dim, ScalarOption option, bool errorTrack>
    class IntegrateSolver<Simpson, dim, option, errorTrack> {

    };

    template<ScalarOption option, bool errorTrack>
    class IntegrateSolver<Simpson, 1, option, errorTrack> {
        Scalar<option, false> stepSize;
    public:
        explicit IntegrateSolver(Scalar<option, false> stepSize);
        /* Operations */
        Scalar<option, errorTrack> solve(const Integrate<1, option, errorTrack>& i) const;
    };
    //////////////////////////////////Tanh_Sinh//////////////////////////////////
    template<size_t dim, ScalarOption option, bool errorTrack>
    class IntegrateSolver<Tanh_Sinh, dim, option, errorTrack> {

    };

    template<ScalarOption option, bool errorTrack>
    class IntegrateSolver<Tanh_Sinh, 1, option, errorTrack> {
        Scalar<option, false> stepSize;
        //Points count on one side of the x-axis. There are same points on both sides.
        size_t pointCount;
    public:
        IntegrateSolver(Scalar<option, false> stepSize, size_t pointCount);
        /* Operations */
        Scalar<option, errorTrack> solve(const Integrate<1, option, errorTrack>& i) const;
    };
}

#include "IntegrateSolverImpl.h"

#endif