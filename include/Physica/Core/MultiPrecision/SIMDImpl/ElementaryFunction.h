/*
 * Copyright 2023 WeiBo He.
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

namespace Physica::Core {
    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> abs(const SIMD<ScalarType, Size>& p) {
        return SIMD<ScalarType, Size>(abs(p.getImpl()));
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> square(const SIMD<ScalarType, Size>& p) {
        return p * p;
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> sqrt(const SIMD<ScalarType, Size>& p) {
        return SIMD<ScalarType, Size>(sqrt(p.getImpl()));
    }

    template<class ScalarType, size_t Size>
    [[nodiscard]] inline SIMD<ScalarType, Size> cbrt(const SIMD<ScalarType, Size>& p) {
        return SIMD<ScalarType, Size>(cbrt(p.getImpl()));
    }

    template<class ScalarType, size_t Size>
    inline static void sincos(const SIMD<ScalarType, Size>& x, SIMD<ScalarType, Size>& s, SIMD<ScalarType, Size>& c) {
        sincos(x.getImpl(), s.getImpl(), c.getImpl());
    }
}
