/*
 * Copyright 2021-2026 Weibo He.
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

namespace Physica {
    /**
     * Reference
     * [1] Test Functions and Datasets; http://www.sfu.ca/~ssurjano/index.html
     */
    template<Scalar T>
    T ackley(const VectorND<T>& x) noexcept {
        return T(20) * (T(1) - exp(T(-0.2) * x.norm2())) - exp(cos(T(2) * MathConst<T>::pi * x).mean()) + exp(T(1));
    }

    template<Scalar T>
    T rosenbrock(const VectorND<T>& x) noexcept {
        assert(x.getLength() >= 2);
        T result = 0;
        for (size_t i = 0; i < x.getLength() - 1; ++i)
            result += fma(T(100), square(fma(x[i], x[i], -x[i + 1])), square(T(1) - x[i]));
        return result;
    }

    template<Scalar T>
    T forrester(T x) noexcept {
        assert(T(0) <= x && x <= T(1) && "[Error]: Out of expected domain");
        return square(fma(T(6), x, T(-2))) * sin(fma(T(12), x, T(-4)));
    }

    template<Scalar T>
    T func1(const VectorND<T>& v) noexcept {
        const T& x = v[0];
        const T& y = v[1];
        const T& z = v[2];
        return x * x + y * y + z * z - x * y - x * z - y * z;
    }

    template<Scalar T>
    T func2(const VectorND<T>& v) noexcept {
        const T& x = v[0];
        const T& y = v[1];
        const T& z = v[2];
        const T term1 = x + y;
        const T term2 = y + z;
        const T term3 = x + z;
        return (reciprocal(square(term1)) + reciprocal(square(term2)) + reciprocal(square(term3))) * T(x * y + x * z + y * z);
    }
}
