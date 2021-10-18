/*
 * Copyright 2021 WeiBo He.
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
    template<ScalarOption option, bool errorTrack>
    Scalar<option, errorTrack> lnGamma(const Scalar<option, errorTrack>& s);
    /**
     * Reference:
     * [1] Hochstrasser, U. W. “Orthogonal Polynomials.” Handbook of Mathematical Functions with
     * Formulas, Graphs, and Mathematical Tables. (M. Abramowitz and I. A. Stegun, eds.). New York: Dover, 1972.789
     */
    template<class ScalarType>
    ScalarType hermiteH(unsigned int n, const ScalarBase<ScalarType>& x) {
        const ScalarType square_x = square(x.getDerived());
        ScalarType a = ScalarType::One();
        ScalarType b = ScalarType::Two();

        const unsigned int k = n / 2;
        ScalarType float_m = ScalarType(k);
        const bool isEven = n % 2 == 0;
        const ScalarType term = isEven ? -ScalarType::One() : ScalarType::One();
        ScalarType c = float_m * (ScalarType::Two() * float_m + term);

        unsigned int m = k;
        while (m > 0) {
            --m;
            float_m -= ScalarType::One();
            a = ScalarType::One() - (b / c) * (square_x * a);
            b += ScalarType::Two();
            c = float_m * (ScalarType::Two() * float_m + term);
        }

        ScalarType result = exp(lnGamma(ScalarType(n + 1)) - lnGamma(ScalarType(k + 1)));
        if (!isEven)
            result *= ScalarType::Two() * x.getDerived();
        result *= (k % 2 == 0) ? a : -a;
        return result;
    }
}
