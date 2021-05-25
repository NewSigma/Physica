/*
 * Copyright 2021 WeiBo He.
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

#include "Physica/SystemBits.h"
#include "Physica/Core/MultiPrecision/MultiPrecisionType.h"

namespace Physica::Core::Internal {
    /**
     * This function is shared by \class Scalar and \class Integer.
     */
    inline double convertDoubleImpl(int length, int power, MPUnit* __restrict byte) {
            double_extract extract{0};
            extract.sign = length < 0;

            const auto size = std::abs(length);
            const auto zeroCount = countLeadingZeros(byte[size - 1]); //Optimize: (size - 1) and (size > 1) is used several times
            //Using long to avoid overflow.
            const long exp = power * __WORDSIZE + static_cast<long>(MPUnitWidth - zeroCount) - 1 + 1023;
            if(exp >= 2047) {
                extract.high = extract.low = 0;
                extract.exp = 2047;
                return extract.value;
            }
            if(exp <= 0) {
                return 0.0;
            }
            extract.exp = exp;

            auto temp = byte[size - 1] << (zeroCount + 1);
        #ifdef PHYSICA_64BIT
            extract.high = temp >> 44U;
            if(zeroCount <= 11) {
                extract.low = temp << 20U >> 32U;
            }
            else {
                if(zeroCount <= 44 - 1) {
                    extract.low = temp << 20U >> 32U;
                    if(size > 1)
                        extract.low += byte[size - 2] >> (64 - (32 - (64 - 20 - zeroCount - 1)));
                }
                else {
                    if(size > 1) {
                        extract.high += byte[size - 2] >> (64 - (20 - (64 - zeroCount - 1)));
                        extract.low = byte[size - 2] << (20 - (64 - zeroCount - 1)) >> 32U;
                    }
                }
            }
        #endif
        #ifdef PHYSICA_32BIT
            extract.high = temp >> 12U;
            if(zeroCount <= 11) {
                extract.low = temp << 20U;
                if(size > 1)
                    extract.low = byte[size - 1] >> (32 - 20 - zeroCount - 1);
            }
            else {
                if(size > 1) {
                    extract.high += byte[size - 1] >> (32 - (zeroCount + 1 - 12));
                    extract.low = byte[size - 1] << (zeroCount + 1 - 12);
                }
                if(size > 2)
                    extract.low += byte[size - 2] >> (32 - (zeroCount + 1 - 12));
            }
        #endif
            return extract.value;
    }
}