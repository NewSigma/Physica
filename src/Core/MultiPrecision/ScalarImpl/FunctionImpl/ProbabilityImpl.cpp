/*
 * Copyright 2023 Weibo He.
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
#include "Physica/Core/MultiPrecision/Scalar.h"

namespace Physica::Core {
    Scalar<MultiPrecision> floor(const Scalar<MultiPrecision>& s) noexcept {
        if(s.isInteger())
            return Scalar(s);
        const auto size = s.getSize();
        const auto power = s.getPower();
        const auto power_1 = power + 1;
        auto length = size > power_1 ? power_1 : size;
        length = s.isNegative() ? -length : length;
        Scalar<MultiPrecision> result(length, power);
        for(int i = 0; i < length; ++i)
            result.setByte(i, s[i]);
        return result;
    }
}
