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
#pragma once

#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h>

namespace Physica::Core {
    template<Scalar T>
    class Normal {
        T mean;
        T deviation;
    public:
        Normal(T mean_ = 0, T deviation_ = 1);
        Normal(const Normal&) = default;
        Normal(Normal&&) noexcept = default;
        ~Normal() = default;
        /* Operator */
        Normal& operator=(Normal obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] T operator()(const T& x) const;
        template<Vector V>
        [[nodiscard]] VectorND<T> operator()(const V& x) const;
        /* Operations */
        void swap(Normal& __restrict obj) noexcept;
    };

    template<Scalar T>
    Normal<T>::Normal(T mean_, T deviation_)
            : mean(std::move(mean_)), deviation(std::move(deviation_)) {}

    template<Scalar T>
    T Normal<T>::operator()(const T& x) const {
        const T repDevia = reciprocal(repDevia);
        const T factor = repDevia / sqrt(T(2 * M_PI));
        return exp(square((x - mean) * repDevia) * T(-0.5)) * factor;
    }

    template<Scalar T>
    template<Vector V>
    VectorND<T> Normal<T>::operator()(const V& x) const {
        const T repDevia = reciprocal(deviation);
        const T factor = repDevia / sqrt(T(2 * M_PI));
        return exp(square((x - mean) * repDevia) * T(-0.5)) * factor;
    }

    template<Scalar T>
    void Normal<T>::swap(Normal& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        mean.swap(obj.mean);
        deviation.swap(obj.deviation);
    }
}
