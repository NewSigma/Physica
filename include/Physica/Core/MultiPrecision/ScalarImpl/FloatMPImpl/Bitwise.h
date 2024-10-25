/*
 * Copyright 2019-2024 Weibo He.
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

#include <Physica/Core/MultiPrecision/MultiPrecisionType.h>

namespace Physica::Core {
    inline unsigned int countLeadingZeros(MPUnit n) noexcept {
        if(n == 0)
            return MPUnitWidth;

        MPUnit count;
        if constexpr (UseASM()) {
            asm volatile (
                "bsrq %1, %0\n\t"
                : "=r" (count)
                : "rm" (n)
            );
            (count) ^= 63U;
        }
        else {
            count = 0;
            while((n & MPUnitHighestBitMask) == 0) {
                ++count;
                n <<= 1U;
            }
        }
        return count;
    }

    inline unsigned int countBackZeros(unsigned long n) noexcept {
        if(n == 0)
            return MPUnitWidth;

        MPUnit count;
        if constexpr (UseASM()) {
            asm volatile (
                "bsfq %1, %0\n\t"
                : "=r" (count)
                : "rm" (n)
            );
        }
        else {
            count = 0;
            while((n & 1U) == 0) {
                ++count;
                n >>= 1U;
            }
        }
        return count;
    }

    inline int countOnes(MPUnit n) noexcept {
        int count=0 ;
        while(n != 0) {
            n &= (n - 1);
            count += 1;
        }
        return count;
    }
}
