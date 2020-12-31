/*
 * Copyright 2020 WeiBo He.
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
#ifndef PHYSICA_PROBABILITYFUNCTION_H
#define PHYSICA_PROBABILITYFUNCTION_H
/*!
 * This file is part of implementations of \Scalar.
 * Do not include this header file, include Scalar.h instead.
 */
namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> randomScalar();
    
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> randomScalar(
            const Scalar<type, errorTrack>& lowerBound,
            const Scalar<type, errorTrack>& upperBound);

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> floor(const Scalar<type, errorTrack>& s);
    
    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> ceil(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> arrangement(const Scalar<type, errorTrack>& s1, const Scalar<type, errorTrack>& s2);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> combination(const Scalar<type, errorTrack>& s1, const Scalar<type, errorTrack>& s2);
}

#include "FunctionImpl/ProbabilityImpl.h"

#endif