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

#include "Physica/Core/MultiPrecision/Scalar.h"
#include "Physica/Core/Math/NumberTheory/NumberTheory.h"

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> lnGamma(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> gamma(const Scalar<type, errorTrack>& s);

    template<ScalarType type, bool errorTrack>
    inline Scalar<type, errorTrack> beta(const Scalar<type, errorTrack>& s1, const Scalar<type, errorTrack>& s2);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> gammaP(const Scalar<type, errorTrack>& a, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> gammaQ(const Scalar<type, errorTrack>& a, const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> erf(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> erfc(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ0(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJ1(const Scalar<type, errorTrack>& x);

    template<ScalarType type, bool errorTrack>
    Scalar<type, errorTrack> besselJn(const Integer& n, const Scalar<type, errorTrack>& x);
}

#include "SpetialFunctionsImpl/Gamma.h"
#include "SpetialFunctionsImpl/Bessel.h"
