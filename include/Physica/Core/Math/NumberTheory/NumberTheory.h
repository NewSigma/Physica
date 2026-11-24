/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <utility>
#include "Physica/Core/Math/Discrete/Combination.h"

namespace Physica {
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
            result -= binomial(i, temp) * bernoulli(temp) / (i - temp + 1);
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
                if (i1 < i2)
                    std::swap(i1, i2);
            }
            i1 = i1 << shift;
        }
        return i1;
    }

    template<class IntType, bool IsBigInt>
    IntType lcm(IntType i1, IntType i2) noexcept {
        IntType result = i1 * i2;
        result /= gcd<IntType, IsBigInt>(std::move(i1), std::move(i2));
        return result;
    }
}
