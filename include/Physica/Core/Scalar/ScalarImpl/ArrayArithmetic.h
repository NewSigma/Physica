/*
 * Copyright 2020-2026 Weibo He.
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

#include <memory>
#include "Physica/Core/Scalar/Scalar.h"

namespace Physica {
    /**
     * \len is the length of \from. Length of \to should not less than \len. length of \result should be \len at least.
     */
    PHYSICA_API MPUnit addArrWithArr(MPUnit* __restrict result, const MPUnit* __restrict from, const MPUnit* __restrict to, size_t len) noexcept;
    /**
     * \len is the length of \from. Length of \to should not less than \len.
     */
    PHYSICA_API MPUnit addArrWithArrEq(const MPUnit* __restrict from, MPUnit* __restrict to, size_t len) noexcept;
    /*
     * Calculate n1 - n2, the result must be one word.
     */
    PHYSICA_API void sub2WordByWord(MPUnit& high, MPUnit& low, MPUnit n1_high, MPUnit n1_low, MPUnit n2) noexcept;
    //Another version of sub2WordByWord()
    PHYSICA_API void sub2WordByWord(MPUnit& high, MPUnit& low, MPUnit n2) noexcept;
    /*
     * Calculate arr1 - arr2. len is the length of arr1 and arr2.
     * If arr1 >= arr2 the function will return 0, if not, the calculation is failed and return true.
     */
    PHYSICA_API bool subArrByArr(MPUnit* __restrict result, const MPUnit* __restrict arr1, const MPUnit* __restrict arr2, size_t len) noexcept;
    //Another version of subArrByArr(), calculate arr1 -= arr2.
    PHYSICA_API bool subArrByArrEq(MPUnit* __restrict arr1, const MPUnit* __restrict arr2, size_t len) noexcept;
    /**
     * This is simplified version of mulWordByWord(), which get the high Unit only.
     * It is slightly faster than mulWordByWord() if we are interested in the high Unit only.
     */
    PHYSICA_API MPUnit mulWordByWordHigh(MPUnit n1, MPUnit n2) noexcept;
    /**
     * On 64 bits machine(similar to 32 bit machine):
     * n1 * n2 = product(16 bytes) = carry(high 8 bytes) + ReturnValue(low bytes)
     */
    PHYSICA_API void mulWordByWord(MPUnit& high, MPUnit& low, MPUnit n1, MPUnit n2) noexcept;
    /**
     * Multiply the array @param arr with @param n. Write the result to array @param result.
     * Length of result should at least as long as arr.
     */
    PHYSICA_API MPUnit mulArrByWord(MPUnit* __restrict result, const MPUnit* __restrict arr, size_t length, MPUnit n) noexcept;
    /**
     * Multiply the array @param arr with @param n. Add the result to array @param result.
     * Length of result should at least as long as arr.
     */
    PHYSICA_API MPUnit mulAddArrByWord(MPUnit* __restrict result, const MPUnit* __restrict arr, size_t length, MPUnit n) noexcept;
    /**
     * Multiply the array @param arr with @param n. Subtract the result from array @param result.
     * Length of result should at least as long as arr.
     */
    PHYSICA_API MPUnit mulSubArrByWord(MPUnit* __restrict result, const MPUnit* __restrict arr, size_t length, MPUnit n) noexcept;
    /*
     * Return the precomputed reciprocal.
     * Reference:
     * [1] T. Granlund and N. M¨oller, “Division of integers large and small”, to appear.
     */
    PHYSICA_API MPUnit getInverse(MPUnit unit) noexcept;
    /*
     * Calculate (high, low) / divisor.
     * Assume high < divisor and divisor >= 2^(MPUnitWidth - 1).
     *
     * A full word here indicates that the highest bit of divisor is set.
     *
     * Reference:
     * [1] T. Granlund and N. M¨oller, “Division of integers large and small”, to appear.
     */
    PHYSICA_API void div2WordByFullWord(MPUnit& quotient, MPUnit& remainder, MPUnit high, MPUnit low, MPUnit divisor) noexcept;
    /*
     * This is a simplified version of div2WordByFullWord(), which returns the quotient only.
     */
    PHYSICA_API MPUnit div2WordByFullWordQ(MPUnit high, MPUnit low, MPUnit divisor) noexcept;
    /*
     * This is a simplified version of div2WordByFullWord(), which returns the remainder only.
     */
    PHYSICA_API MPUnit div2WordByFullWordR(MPUnit high, MPUnit low, MPUnit divisor) noexcept;
    /*
     * Requirement:
     * The quotient must be one word.
     * len is the length of divisor, length of dividend should equals to len + 1.
     * FullArr here indicates that the highest bit of divisor is set.
     * dividend[len] < divisor[len - 1].
     * len >= 1 to avoid invalid read.
     *
     * Reference:
     * [1] MaTHmu Project Group.计算机代数系统的数学原理[M].Beijing: TsingHua University Press, 2009:4-8
     */
    PHYSICA_API MPUnit divArrByFullArrWith1Word(const MPUnit* __restrict dividend, const MPUnit* __restrict divisor, size_t len) noexcept;
    // operator<<
    [[nodiscard]] PHYSICA_API std::unique_ptr<MPUnit[]> byteLeftShift(const MPUnit* __restrict byte, unsigned int length, unsigned int shift) noexcept;
    // operator>>
    [[nodiscard]] PHYSICA_API std::unique_ptr<MPUnit[]> byteRightShift(const MPUnit* __restrict byte, size_t length, size_t shift) noexcept;
    // operator<<=
    PHYSICA_API void byteLeftShiftEq(MPUnit* __restrict byte, unsigned int length, unsigned int shift) noexcept;
    // operator>>=
    PHYSICA_API void byteRightShiftEq(MPUnit* __restrict byte, size_t length, size_t shift) noexcept;
}
