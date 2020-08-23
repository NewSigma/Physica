/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_BITWISE_H
#define PHYSICA_BITWISE_H

#include "Physica/SystemBits.h"

namespace Physica::Core {
    //Possibly use asm to speed up.
    inline unsigned int countLeadingZeros(ScalarUnit n) {
        if(n == 0)
            return ScalarUnitWidth;
        ScalarUnit count;
    #if UseASM
        asm volatile (
                "bsrq %1, %0\n\t"
                : "=r" (count)
                : "rm" (n)
        );
        (count) ^= 63U;
    #else
        count = 0;
        while((n & ScalarUnitHighestBitMask) == 0) {
            ++count;
            n <<= 1U;
        }
    #endif
        return count;
    }

    inline unsigned int countBackZeros(unsigned long n) {
        if(n == 0)
            return ScalarUnitWidth;
        ScalarUnit count;
    #if UseASM
        asm volatile (
        "bsfq %1, %0\n\t"
        : "=r" (count)
        : "rm" (n)
        );
    #else
        count = 0;
        while((n & 1U) == 0) {
            ++count;
            n >>= 1U;
        }
    #endif
        return count;
    }
}

#endif