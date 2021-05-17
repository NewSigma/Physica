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

#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core {
    template<class T, size_t Length, size_t Capacity>
    std::pair<T, T> linearFit(const Utils::Array<T, Length, Capacity>& x, const Utils::Array<T, Length, Capacity>& y) {
        assert(x.getLength() == y.getLength());
        const size_t length = x.getLength();
        T xy_sum(0);
        T x_square_sum(0);
        T x_sum(0);
        T y_sum(0);
        for (size_t i = 0; i < length; ++i) {
            const auto& x_i = x[i];
            const auto& y_i = y[i];
            xy_sum += x_i * y_i;
            y_sum += y_i;
            x_square_sum += square(x_i);
            x_sum += x_i;
        }
        const T length_1 = reciprocal(T(length));
        T numerator = xy_sum - x_sum * y_sum * length_1;
        T denominator = x_square_sum - square(x_sum) * length_1;
        T k = numerator / denominator;
        return {k, (y_sum - k * x_sum) * length_1};
    }
}