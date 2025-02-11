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
        template<Vector V>
        [[nodiscard]] VectorND<T> operator()(const V& x) const;
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
    template<Vector V>
    VectorND<T> LogNormal<T>::operator()(const V& x) const {
        return divide(normal(ln(x)), x);
    }

    template<Scalar T>
    void LogNormal<T>::swap(LogNormal& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        normal.swap(obj.normal);
    }
}
