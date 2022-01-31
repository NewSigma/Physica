/*
 * Copyright 2021-2022 WeiBo He.
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

    template<class ScalarType>
    ScalarType hermiteH(unsigned int n, const ScalarBase<ScalarType>& x) {
        using std::swap;
        if (n == 0)
            return ScalarType::One();
        const ScalarType double_x = ScalarType::Two() * x.getDerived();
        if (n == 1)
            return double_x;

        ScalarType old_H = ScalarType::One();
        ScalarType H = double_x;
        ScalarType float_i = ScalarType(1);
        for (unsigned int i = 1; i != n; ++i) {
            old_H = double_x * H - float_i * old_H * ScalarType::Two();
            swap(old_H, H);
            float_i += ScalarType::One();
        }
        return H;
    }
}
