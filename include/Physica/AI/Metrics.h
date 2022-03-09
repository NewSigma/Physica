/*
 * Copyright 2022 WeiBo He.
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

#include <climits>

namespace Physica::AI {
    /**
     * \returns polarization rate, the lower the better, minus value means overfitting
     */
    template<class ScalarType>
    ScalarType polarization(ScalarType train_loss, ScalarType valid_loss) {
        const ScalarType total = train_loss + valid_loss;
        const ScalarType delta = train_loss - valid_loss;
        if (abs(delta) >= total * ScalarType(std::numeric_limits<ScalarType>::epsilon())):
            return ScalarType(0);
        return delta / total * ScalarType(2);
    }
}