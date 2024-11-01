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
    inline int countLeadingZeros(MPUnit n) noexcept {
        if(n == 0)
            return MPUnitWidth;

        MPUnit count;
        if constexpr (UseASM() && !IsMSVC()) {
        #ifdef __GNUC__
            asm volatile (
                "bsrq %1, %0\n\t"
                : "=r" (count)
                : "rm" (n)
            );
            (count) ^= 63U;
        #else
            noImpl();
        #endif
        }
        else {
            count = 0;
            while((n & MPUnitHighestBitMask) == 0) {
                ++count;
                n <<= 1U;
            }
        }
        return int(count);
    }

    inline int countBackZeros(unsigned long n) noexcept {
        if(n == 0)
            return MPUnitWidth;

        MPUnit count;
        if constexpr (UseASM() && !IsMSVC()) {
        #ifdef __GNUC__
            asm volatile (
                "bsfq %1, %0\n\t"
                : "=r" (count)
                : "rm" (n)
            );
        #else
            noImpl();
        #endif
        }
        else {
            count = 0;
            while((n & 1U) == 0) {
                ++count;
                n >>= 1U;
            }
        }
        return int(count);
    }

    inline int countOnes(MPUnit n) noexcept {
        int count=0 ;
        while(n != 0) {
            n &= (n - 1);
            count += 1;
        }
        return count;
    }
    /**
     * This function is shared by \class Scalar and \class Integer.
     */
    inline double convertDoubleImpl(int length, int power, MPUnit* __restrict byte) {
            double_extract extract{0};
            extract.sign = length < 0;

            const auto size = std::abs(length);
            const auto zeroCount = countLeadingZeros(byte[size - 1]); //Optimize: (size - 1) and (size > 1) is used several times
            //Using long to avoid overflow.
            const long exp = power * int(MPUnitWidth) + static_cast<long>(MPUnitWidth - zeroCount) - 1 + 1023;
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
            if constexpr (PhysicaWordSize == 64) {
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
            }
            else {
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
            }
            return extract.value;
    }
}
