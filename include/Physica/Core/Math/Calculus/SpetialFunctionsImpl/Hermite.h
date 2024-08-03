/*
 * Copyright 2021-2023 Weibo He.
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
    template<ScalarOption Option, bool errorTrack>
    Scalar<Option> lnGamma(const Scalar<Option>& s);

    template<class ScalarType>
    ScalarType hermiteH(unsigned int n, const ScalarBase<ScalarType>& x) {
        using std::swap;
        if (n == 0)
            return ScalarType(1);
        const ScalarType double_x = ScalarType(2) * x.getDerived();
        if (n == 1)
            return double_x;

        ScalarType old_H = ScalarType(1);
        ScalarType H = double_x;
        ScalarType float_i = ScalarType(1);
        for (unsigned int i = 1; i != n; ++i) {
            old_H = double_x * H - float_i * old_H * ScalarType(2);
            swap(old_H, H);
            float_i += ScalarType(1);
        }
        return H;
    }
}
