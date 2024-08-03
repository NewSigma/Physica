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

namespace Physica::Core {
    template<class ScalarType>
    class LogNormal {
        Normal<ScalarType> normal;
    public:
        LogNormal(ScalarType logMean, ScalarType logDevia);
        LogNormal(const LogNormal&) = default;
        LogNormal(LogNormal&&) noexcept = default;
        ~LogNormal() = default;
        /* Operators */
        LogNormal& operator=(LogNormal obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] ScalarType operator()(const ScalarType& x) const;
        template<class VectorType>
        [[nodiscard]] Vector<ScalarType> operator()(const RValueVector<VectorType>& x) const;
        /* Operations */
        void swap(LogNormal& __restrict obj) noexcept;
    };

    template<class ScalarType>
    LogNormal<ScalarType>::LogNormal(ScalarType logMean, ScalarType logDevia)
            : normal(std::move(logMean), std::move(logDevia)) {}

    template<class ScalarType>
    ScalarType LogNormal<ScalarType>::operator()(const ScalarType& x) const {
        return normal(ln(x)) / x;
    }

    template<class ScalarType>
    template<class VectorType>
    Vector<ScalarType> LogNormal<ScalarType>::operator()(const RValueVector<VectorType>& x) const {
        return divide(normal(ln(x)), x);
    }

    template<class ScalarType>
    void LogNormal<ScalarType>::swap(LogNormal& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        normal.swap(obj.normal);
    }
}
