/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

namespace Physica {
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
        [[nodiscard]] VectorND<T> operator()(const Vector auto& x) const;
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
    VectorND<T> Normal<T>::operator()(const Vector auto& x) const {
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
