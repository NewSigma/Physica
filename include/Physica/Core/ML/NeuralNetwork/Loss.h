/*
 * Copyright 2023-2025 Weibo He.
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

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica::Core {
    /**
     * \returns polarization rate, the lower the better, minus value means overfitting
     */
    template<Scalar T>
    T polar_rate(T train_loss, T valid_loss) {
        const T total = train_loss + valid_loss;
        const T delta = train_loss - valid_loss;
        if (abs(delta) * T(std::numeric_limits<T>::epsilon()) >= total)
            return T(0);
        return delta / total * T(2);
    }

    template<Scalar T>
    T mixed_loss(T train_loss, T valid_loss) {
        return std::max(train_loss, valid_loss) + abs(train_loss - valid_loss);
    }
}
