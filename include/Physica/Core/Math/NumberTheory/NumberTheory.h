/*
 * Copyright 2020-2024 WeiBo He.
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

#include "Physica/Core/Math/Discrete/Combination.h"

namespace Physica::Core {
    class Integer;
    /**
     * Optimize: Make use of DP and the fact that $C_n^m = C_n^(n - m)$
     */
    template<class IntType>
    IntType bernoulli(const IntType& i) {
        IntType result = 1;
        if(i.isZero())
            return result;
        IntType temp = 0;
        for(; temp < i; ++temp)
            result -= combination(i, temp) * bernoulli(temp) / (i - temp + 1);
        return result;
    }
    /**
     * \tparam IsBigInt
     * Use euclidean algorithm if true, use decreases technique if false.
     * 
     * Decreases technique runs faster because subtract is faster than division
     */
    template<class IntType, bool IsBigInt>
    IntType gcd(IntType i1, IntType i2) {
        if (i1 < i2)
            std::swap(i1, i2);
        if constexpr (IsBigInt) {
            while (!i2.isZero()) {
                i1 = i1 % i2;
                std::swap(i1, i2);
            }
        }
        else {
            int shift = 0;
            constexpr bool IsMultiInt = std::is_same<IntType, Integer>::value;
            if constexpr (IsMultiInt) {
                while(i1.isEven() && i2.isEven()) {
                    assert(shift < INT_MAX && "[Error]: Params are too large, use euclidean algorithm instead");
                    i1 >>= 1;
                    i2 >>= 1;
                    ++shift;
                }
            }
            else {
                while((i1 % IntType(2) == 0) && (i2 % IntType(2) == 0)) {
                    assert(shift < INT_MAX && "[Error]: Params are too large, use euclidean algorithm instead");
                    i1 >>= 1;
                    i2 >>= 1;
                    ++shift;
                }
            }


            while(i2 != IntType(0)) {
                i1 = i1 - i2;
                std::swap(i1, i2);
            }
            i1 = i1 << shift;
        }
        return i1;
    }

    template<class IntType, bool IsBigInt>
    inline IntType lcm(IntType i1, IntType i2) noexcept {
        IntType result = i1 * i2;
        result /= gcd<IntType, IsBigInt>(std::move(i1), std::move(i2));
        return result;
    }
}
