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

namespace Physica {
    template<ScalarOption Option, bool errorTrack>
    Real<Option> lnGamma(const Real<Option>& s);

    template<Scalar T>
    T hermiteH(unsigned int n, const T& x) {
        if (n == 0)
            return T(1);
        const T double_x = T(2) * x;
        if (n == 1)
            return double_x;

        T old_H = T(1);
        T H = double_x;
        T float_i = T(1);
        for (unsigned int i = 1; i != n; ++i) {
            old_H = double_x * H - float_i * old_H * T(2);
            old_H.swap(H);
            float_i += T(1);
        }
        return H;
    }
}
