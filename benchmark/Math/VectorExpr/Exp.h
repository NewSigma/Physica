/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
    template<Scalar T>
    VectorND<T> makeData(size_t size) {
        using RandomSource = Random<MCG, std::mt19937::default_seed>;
        const Array<float64, 8> factors{1, 2, 10, 100, -1, -2, -10, -100}; // Covers wide numeric range

        RandomSource::getInstance().reseed();
        VectorND<T> data = VectorND<T>::template random_uniform<RandomSource>(size);
        for (size_t i = 0; i < size; ++i)
            data[i] *= factors[i % factors.size()];
        std::ranges::shuffle(data, RandomSource::getInstance());
        return data;
    }
}
