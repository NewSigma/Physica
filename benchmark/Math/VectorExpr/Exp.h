/*
 * Copyright 2025 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    template<Scalar T>
    VectorND<T> makeData(size_t size) {
        using RandomSource = Random<>;
        const Array<float64, 8> factors{1, 2, 10, 100, -1, -2, -10, -100}; // Covers wide numeric range

        VectorND<T> data = VectorND<T>::template random_uniform<RandomSource>(size);
        for (size_t i = 0; i < size; ++i)
            data[i] *= factors[i % factors.size()];
        std::ranges::shuffle(data, RandomSource::getInstance());
        return data;
    }
}
