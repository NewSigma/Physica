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
    inline Scalar<option> floor(const Scalar<option>& s);
    
    template<ScalarOption option>
    inline Scalar<option> ceil(const Scalar<option>& s);

    template<ScalarOption option>
    Scalar<option> arrangement(const Scalar<option>& s1, const Scalar<option>& s2);

    template<ScalarOption option>
    Scalar<option> combination(const Scalar<option>& s1, const Scalar<option>& s2);
}

#include "FunctionImpl/ProbabilityImpl.h"
