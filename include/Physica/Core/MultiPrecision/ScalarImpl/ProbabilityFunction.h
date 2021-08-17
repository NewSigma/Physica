/*
 * Copyright 2020-2021 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType>
    inline ScalarType randomScalar();
    
    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> randomScalar(
            const Scalar<option, errorTrack>& lowerBound,
            const Scalar<option, errorTrack>& upperBound);

    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> floor(const Scalar<option, errorTrack>& s);
    
    template<ScalarOption option, bool errorTrack>
    inline Scalar<option, errorTrack> ceil(const Scalar<option, errorTrack>& s);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> arrangement(const Scalar<option, errorTrack>& s1, const Scalar<option, errorTrack>& s2);

    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> combination(const Scalar<option, errorTrack>& s1, const Scalar<option, errorTrack>& s2);
}

#include "FunctionImpl/ProbabilityImpl.h"
