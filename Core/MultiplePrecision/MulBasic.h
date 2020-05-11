/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_MULBASIC_H
#define PHYSICA_MULBASIC_H

#include "Core/Header/SystemBits.h"

namespace Physica::Core {
    extern const NumericalUnit numericalUnitLowMask;
    /*
     * This is simplified version of mulWordByWord(), which get the high Unit only.
     * It is slightly faster than mulWordByWord() if we are interested in the high Unit only.
     *
     * Possibly use asm to speed up. Depending on the platform.
     */
    inline NumericalUnit mulWordByWordHigh(NumericalUnit n1, NumericalUnit n2) {
        unsigned long n1_low = n1 & numericalUnitLowMask;
        unsigned long n1_high = n1 >> (64U / 2U);
        unsigned long n2_low = n2 & numericalUnitLowMask;
        unsigned long n2_high = n2 >> (64U / 2U);

        auto ll = n1_low * n2_low;
        auto lh = n1_low * n2_high;
        auto hl = n1_high * n2_low;
        auto hh = n1_high * n2_high;

        lh += ll >> (64U / 2U);
        lh += hl;
        hh += static_cast<unsigned long>(lh < hl) << (64U / 2U);

        return hh + (lh >> (64U / 2U));
    }
/*
 * On 64 bits machine(similar to 32 bit machine):
 * n1 * n2 = product(16 bytes) = carry(high 8 bytes) + ReturnValue(low bytes)
 *
 * Possibly use asm to speed up. Depending on the platform.
 */
    inline void mulWordByWord(NumericalUnit& high, NumericalUnit& low, NumericalUnit n1, NumericalUnit n2) {
        NumericalUnit n1_low = n1 & numericalUnitLowMask;
        NumericalUnit n1_high = n1 >> (NumericalUnitWidth / 2U);
        NumericalUnit n2_low = n2 & numericalUnitLowMask;
        NumericalUnit n2_high = n2 >> (NumericalUnitWidth / 2U);

        auto ll = n1_low * n2_low;
        auto lh = n1_low * n2_high;
        auto hl = n1_high * n2_low;
        auto hh = n1_high * n2_high;

        lh += ll >> (NumericalUnitWidth / 2U);
        lh += hl;
        hh += static_cast<NumericalUnit>(lh < hl) << (NumericalUnitWidth / 2U);
        high = hh + (lh >> (NumericalUnitWidth / 2U));
        low = (lh << (NumericalUnitWidth / 2U)) + (ll & numericalUnitLowMask);
    }
}

#endif