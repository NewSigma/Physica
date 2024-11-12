/*
 * Copyright 2024 Weibo He.
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
    /**
     * References:
     * [1] add-on; https://github.com/vectorclass/add-on
     */
    template<class T, size_t Size>
    inline SIMD<Complex<T>, Size> sqrt(const SIMD<Complex<T>, Size>& x) {
        using RealPack = SIMD<T, Size * 2>;
        const RealPack x1 = x.asReal();
        const RealPack t1 = x1 * x1;
        RealPack t2;
        if constexpr (T::Option == Float32)
            t2 = t1.template shuffle<1, 0, 3, 2>();
        else {
            static_assert(T::Option == Float64, "[Error]: Not implemented");
            if constexpr (Size == 1)
                t2 = t1.template shuffle<1, 0>();
            else if constexpr (Size == 2)
                t2 = t1.template shuffle<1, 0, 1, 0>();
            else {
                static_assert(Size == 4, "[Error]: Unexpected size");
                t2 = t1.template shuffle<1, 0, 1, 0, 1, 0, 1, 0>();
            }
        }

        const RealPack t3 = sqrt(t1 + t2);
        RealPack t4;
        if constexpr (T::Option == Float32)
            t4 = x1.template shuffle<0, 0, 2, 2>();
        else {
            static_assert(T::Option == Float64, "[Error]: Not implemented");
            if constexpr (Size == 1)
                t4 = x1.template shuffle<0, 0>();
            else if constexpr (Size == 2)
                t4 = x1.template shuffle<0, 0, 0, 0>();
            else {
                static_assert(Size == 4, "[Error]: Unexpected size");
                t4 = x1.template shuffle<0, 0, 0, 0, 0, 0, 0, 0>();
            }
        }

        RealPack signbit;
        if constexpr (Size == 1)
            signbit = RealPack::template makeSignBits<0, 1>();
        else if constexpr (Size == 2)
            signbit = RealPack::template makeSignBits<0, 1, 0, 1>();
        else if constexpr (Size == 4)
            signbit = RealPack::template makeSignBits<0, 1, 0, 1, 0, 1, 0, 1>();
        else {
            static_assert(Size == 8, "[Error]: Unexpected size");
            signbit = RealPack::template makeSignBits<0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1>();
        }
        const RealPack t5 = t3 + (t4 ^ signbit);
        const RealPack t6 = sqrt(t5 * T(0.5));
        const RealPack result = t6 ^ (x1 & signbit);
        return SIMD<Complex<T>, Size>::asComplex(result);
    }
}
