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

    template<IntegrateMethod method, size_t dim, ScalarType type = MultiPrecision, bool errorTrack = true>
    class IntegrateSolver;
    //////////////////////////////////Rectangular//////////////////////////////////
    template<size_t dim, ScalarType type, bool errorTrack>
    class IntegrateSolver<Rectangular, dim, type, errorTrack> {

    };

    template<ScalarType type, bool errorTrack>
    class IntegrateSolver<Rectangular, 1, type, errorTrack> {
        Scalar<type, false> stepSize;
    public:
        explicit IntegrateSolver(Scalar<type, false> stepSize);
        /* Operations */
        Scalar<type, errorTrack> solve(const Integrate<1, type, errorTrack>& i) const;
    };
    //////////////////////////////////Ladder//////////////////////////////////
    template<size_t dim, ScalarType type, bool errorTrack>
    class IntegrateSolver<Ladder, dim, type, errorTrack> {

    };

    template<ScalarType type, bool errorTrack>
    class IntegrateSolver<Ladder, 1, type, errorTrack> {
        Scalar<type, false> stepSize;
    public:
        explicit IntegrateSolver(Scalar<type, false> stepSize);
        /* Operations */
        Scalar<type, errorTrack> solve(const Integrate<1, type, errorTrack>& i) const;
    };
    //////////////////////////////////Simpson//////////////////////////////////
    template<size_t dim, ScalarType type, bool errorTrack>
    class IntegrateSolver<Simpson, dim, type, errorTrack> {

    };

    template<ScalarType type, bool errorTrack>
    class IntegrateSolver<Simpson, 1, type, errorTrack> {
        Scalar<type, false> stepSize;
    public:
        explicit IntegrateSolver(Scalar<type, false> stepSize);
        /* Operations */
        Scalar<type, errorTrack> solve(const Integrate<1, type, errorTrack>& i) const;
    };
    //////////////////////////////////Tanh_Sinh//////////////////////////////////
    template<size_t dim, ScalarType type, bool errorTrack>
    class IntegrateSolver<Tanh_Sinh, dim, type, errorTrack> {

    };

    template<ScalarType type, bool errorTrack>
    class IntegrateSolver<Tanh_Sinh, 1, type, errorTrack> {
        Scalar<type, false> stepSize;
        //Points count on one side of the x-axis. There are same points on both sides.
        size_t pointCount;
    public:
        IntegrateSolver(Scalar<type, false> stepSize, size_t pointCount);
        /* Operations */
        Scalar<type, errorTrack> solve(const Integrate<1, type, errorTrack>& i) const;
    };
}

#include "IntegrateSolverImpl.h"

#endif