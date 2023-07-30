/*
 * Copyright 2020-2023 WeiBo He.
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
    template<ScalarOption option>
    Scalar<option> floor(const Scalar<option>& s) {
        return Scalar<option>(std::floor(s.getTrivial()));
    }

    template<>
    Scalar<MultiPrecision> floor(const Scalar<MultiPrecision>& s);

    template<ScalarOption option>
    inline Scalar<option> ceil(const Scalar<option>& s) {
        return Scalar<option>(std::ceil(s.getTrivial()));
    }
}
