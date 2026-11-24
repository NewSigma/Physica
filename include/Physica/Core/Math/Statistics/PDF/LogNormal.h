/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Normal.h"

namespace Physica {
    template<Scalar T>
    class LogNormal {
        Normal<T> normal;
    public:
        LogNormal(T logMean, T logDevia);
        LogNormal(const LogNormal&) = default;
        LogNormal(LogNormal&&) noexcept = default;
        ~LogNormal() = default;
        /* Operators */
        LogNormal& operator=(LogNormal obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] T operator()(const T& x) const;
        [[nodiscard]] VectorND<T> operator()(const Vector auto& x) const;
        /* Operations */
        void swap(LogNormal& __restrict obj) noexcept;
    };

    template<Scalar T>
    LogNormal<T>::LogNormal(T logMean, T logDevia)
            : normal(std::move(logMean), std::move(logDevia)) {}

    template<Scalar T>
    T LogNormal<T>::operator()(const T& x) const {
        return normal(ln(x)) / x;
    }

    template<Scalar T>
    VectorND<T> LogNormal<T>::operator()(const Vector auto& x) const {
        return divide(normal(ln(x)), x);
    }

    template<Scalar T>
    void LogNormal<T>::swap(LogNormal& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        normal.swap(obj.normal);
    }
}
