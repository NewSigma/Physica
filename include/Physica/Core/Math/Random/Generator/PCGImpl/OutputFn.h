/*
 * PCG Random Number Generation for C++
 *
 * Copyright 2014-2022 Melissa O'Neill <oneill@pcg-random.org>,
 *                     and the PCG Project contributors.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 * Licensed under the Apache License, Version 2.0 (provided in
 * LICENSE-APACHE.txt and at http://www.apache.org/licenses/LICENSE-2.0)
 * or under the MIT license (provided in LICENSE-MIT.txt and at
 * http://opensource.org/licenses/MIT), at your option. This file may not
 * be copied, modified, or distributed except according to those terms.
 *
 * Distributed on an "AS IS" BASIS, WITHOUT WARRANTY OF ANY KIND, either
 * express or implied.  See your chosen license for details.
 *
 * For additional information about the PCG random number generation scheme,
 * visit http://www.pcg-random.org/.
 */

/*
 * Copyright 2025 Weibo He.
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

#include <concepts>
#include <cstdint>
#include "Physica/Core/Utils/Builtin.h"

namespace Physica::Internal {
    enum OutputFnType : char {
        XSH,
        XSH_RS,
        XSH_RR,
        RXS,
        RXS_M,
        RXS_M_XS,
        XSL,
        XSL_RR,
        XSL_RR_RR,
        DXSM,
    };

    template<class OutI, class WorkI, OutputFnType Type>
    class OutputFn {
        using HalfI1 = std::conditional<std::same_as<WorkI, __uint128_t>, uint64_t, void>::type;
        using HalfI2 = std::conditional<std::same_as<WorkI, uint64_t>, uint32_t, HalfI1>::type;
        using HalfI3 = std::conditional<std::same_as<WorkI, uint32_t>, uint16_t, HalfI2>::type;
        using HalfI = std::conditional<std::same_as<WorkI, uint16_t>, uint8_t, HalfI3>::type;

        constexpr static uint8_t WorkBits = sizeof(WorkI) * 8;
        constexpr static uint8_t OutBits = sizeof(OutI) * 8;
        constexpr static uint8_t HalfBits = sizeof(HalfI) * 8;
        constexpr static uint8_t SpareBits = WorkBits - OutBits;
        static_assert(OutBits <= WorkBits / 2, "[Error]: Working int must larger than output int");
    public:
        /* Operators */
        [[nodiscard]] static OutI operator()(WorkI internal) noexcept;
        /* Static members */
        [[nodiscard]] consteval static WorkI multiplier() noexcept;
        [[nodiscard]] consteval static WorkI increment() noexcept;
    private:
        [[nodiscard]] consteval static WorkI multiplierMCG() noexcept;
        [[nodiscard]] consteval static WorkI unmultiplierMCG() noexcept;

        [[nodiscard]] static OutI xsh(WorkI internal) noexcept;
        [[nodiscard]] static OutI xsh_rs(WorkI internal) noexcept;
        [[nodiscard]] static OutI xsh_rr(WorkI internal) noexcept;
        [[nodiscard]] static OutI rxs(WorkI internal) noexcept;
        [[nodiscard]] static OutI rxs_m(WorkI internal) noexcept;
        [[nodiscard]] static OutI rxs_m_xs(WorkI internal) noexcept;
        [[nodiscard]] static OutI xsl(WorkI internal) noexcept;
        [[nodiscard]] static OutI xsl_rr(WorkI internal) noexcept;
        [[nodiscard]] static OutI xsl_rr_rr(WorkI internal) noexcept;
        [[nodiscard]] static OutI dxsm(WorkI internal) noexcept;

        [[nodiscard]] static WorkI rxs_m_xs_inv(WorkI internal) noexcept;
    };

    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::operator()(WorkI internal) noexcept -> OutI {
        switch (Type) {
        case XSH:
            return xsh(internal);
        case XSH_RS:
            return xsh_rs(internal);
        case XSH_RR:
            return xsh_rr(internal);
        case RXS:
            return rxs(internal);
        case RXS_M:
            return rxs_m(internal);
        case RXS_M_XS:
            return rxs_m_xs(internal);
        case XSL:
            return xsl(internal);
        case XSL_RR:
            return xsl_rr(internal);
        case XSL_RR_RR:
            return xsl_rr_rr(internal);
        case DXSM:
            return dxsm(internal);
        default:
            unreachable();
        }
    }

    template<class OutI, class WorkI, OutputFnType Type>
    consteval WorkI OutputFn<OutI, WorkI, Type>::multiplier() noexcept {
        if constexpr (std::same_as<uint8_t, WorkI>)
            return 141U;
        if constexpr (std::same_as<uint16_t, WorkI>)
            return 12829U;
        if constexpr (std::same_as<uint32_t, WorkI>)
            return 747796405U;
        if constexpr (std::same_as<uint64_t, WorkI>)
            return 6364136223846793005ULL;
        else {
            static_assert(std::same_as<__uint128_t, WorkI>);
            return (__uint128_t(2549297995355413924ULL) << 64U) + 4865540595714422341ULL;
        }
    }

    template<class OutI, class WorkI, OutputFnType Type>
    consteval WorkI OutputFn<OutI, WorkI, Type>::increment() noexcept {
        if constexpr (std::same_as<uint8_t, WorkI>)
            return 77U;
        if constexpr (std::same_as<uint16_t, WorkI>)
            return 47989U;
        if constexpr (std::same_as<uint32_t, WorkI>)
            return 2891336453U;
        if constexpr (std::same_as<uint64_t, WorkI>)
            return 1442695040888963407ULL;
        else {
            static_assert(std::same_as<__uint128_t, WorkI>);
            return (__uint128_t(6364136223846793005ULL) << 64U) + 1442695040888963407ULL;
        }
    }

    template<class OutI, class WorkI, OutputFnType Type>
    consteval WorkI OutputFn<OutI, WorkI, Type>::multiplierMCG() noexcept {
        if constexpr (std::same_as<uint8_t, WorkI>)
            return 217U;
        if constexpr (std::same_as<uint16_t, WorkI>)
            return 62169U;
        if constexpr (std::same_as<uint32_t, WorkI>)
            return 277803737U;
        if constexpr (std::same_as<uint64_t, WorkI>)
            return 12605985483714917081ULL;
        else {
            static_assert(std::same_as<__uint128_t, WorkI>);
            return (__uint128_t(17766728186571221404ULL) << 64U) + 12605985483714917081ULL;
        }
    }

    template<class OutI, class WorkI, OutputFnType Type>
    consteval WorkI OutputFn<OutI, WorkI, Type>::unmultiplierMCG() noexcept {
        if constexpr (std::same_as<uint8_t, WorkI>)
            return 105U;
        if constexpr (std::same_as<uint16_t, WorkI>)
            return 28009U;
        if constexpr (std::same_as<uint32_t, WorkI>)
            return 2897767785U;
        if constexpr (std::same_as<uint64_t, WorkI>)
            return 15009553638781119849ULL;
        else {
            static_assert(std::same_as<__uint128_t, WorkI>);
            return (__uint128_t(14422606686972528997ULL) << 64U) + 15009553638781119849ULL;
        }
    }
    /*
     * XSH -- fixed xorshift (to high bits)
     *
     * You shouldn't use this at 64-bits or less.
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::xsh(WorkI internal) noexcept -> OutI {
        constexpr uint8_t TopSpare = 0;
        constexpr uint8_t BottomSpare = SpareBits - TopSpare;
        constexpr uint8_t xshift = (TopSpare + OutBits) / 2;

        internal ^= internal >> xshift;
        OutI result = internal >> BottomSpare;
        return result;
    }
    /**
     * High xorshift, followed by a random shift
     *
     * Fast. A good performer.
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::xsh_rs(WorkI internal) noexcept -> OutI {
        constexpr uint8_t Opbits = [=]() noexcept {
            if (SpareBits >= 64 + 5)
                return 5;
            if (SpareBits >= 32 + 4)
                return 4;
            if (SpareBits >= 16 + 3)
                return 3;
            if (SpareBits >= 4 + 2)
                return 2;
            return SpareBits >= (1 + 1) ? 1 : 0;
        }();
        constexpr uint8_t Mask = (1 << Opbits) - 1;
        constexpr uint8_t MaxRandShift = Mask;
        constexpr uint8_t TopSpare = Opbits;
        constexpr uint8_t BottomSpare = SpareBits - TopSpare;
        constexpr uint8_t xshift = TopSpare + (OutBits + MaxRandShift) / 2;
        uint8_t rshift = Opbits ? uint8_t(internal >> (WorkBits - Opbits)) & Mask : 0;
        internal ^= internal >> xshift;
        return internal >> (BottomSpare - MaxRandShift + rshift);
    }
    /**
     * High xorshift, followed by a random rotate; Fast.
     *
     * A good performer. Slightly better statistically than XSH RS.
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::xsh_rr(WorkI internal) noexcept -> OutI {
        constexpr uint8_t WantedOpbits = [=]() noexcept {
            if (OutBits >= 128)
                return 7;
            if (OutBits >= 64)
                return 6;
            if (OutBits >= 32)
                return 5;
            if (OutBits >= 16)
                return 4;
            return 3;
        }();
        constexpr uint8_t Opbits = SpareBits >= WantedOpbits ? WantedOpbits : SpareBits;
        constexpr uint8_t amplifier = WantedOpbits - Opbits;
        constexpr uint8_t Mask = (1 << Opbits) - 1;
        constexpr uint8_t TopSpare = Opbits;
        constexpr uint8_t BottomSpare = SpareBits - TopSpare;
        constexpr uint8_t xshift = (TopSpare + OutBits) / 2;
        uint8_t rot = Opbits ? uint8_t(internal >> (WorkBits - Opbits)) & Mask : 0;
        uint8_t amprot = (rot << amplifier) & Mask;
        internal ^= internal >> xshift;
        OutI result = OutI(internal >> BottomSpare);
        result = std::rotr(result, amprot);
        return result;
    }
    /**
     * Random xorshift
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::rxs(WorkI internal) noexcept -> OutI {
        constexpr uint8_t Shift = WorkBits - OutBits;
        constexpr uint8_t ExtraShift = (OutBits - Shift) / 2;
        uint8_t rshift = [=]() noexcept {
            if (Shift > (64 + 8))
                return (internal >> (WorkBits - 6)) & 63;
            if (Shift > (32 + 4))
                return (internal >> (WorkBits - 5)) & 31;
            if (Shift > (16 + 2))
                return (internal >> (WorkBits - 4)) & 15;
            if (Shift > (8 + 1))
                return (internal >> (WorkBits - 3)) & 7;
            if (Shift > (4 + 1))
                return (internal >> (WorkBits - 2)) & 3;
            return (Shift > (2 + 1)) ? (internal >> (WorkBits - 1)) & 1 : 0;
        }();

        internal ^= internal >> (Shift + ExtraShift - rshift);
        OutI result = internal >> rshift;
        return result;
    }
    /**
     * Random xorshift, mcg multiply
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::rxs_m(WorkI internal) noexcept -> OutI {
        constexpr uint8_t Opbits = [=]() noexcept {
            if (OutBits >= 128)
                return 6;
            if (OutBits >= 64)
                return 5;
            if (OutBits >= 32)
                return 4;
            if (OutBits >= 16)
                return 3;
            return 2;
        }();
        constexpr uint8_t Shift = WorkBits - OutBits;
        constexpr uint8_t Mask = (1 << Opbits) - 1;
        uint8_t rshift = Opbits ? (internal >> (WorkBits - Opbits)) & Mask : 0;
        internal ^= internal >> (Opbits + rshift);
        internal *= multiplierMCG();
        OutI result = internal >> Shift;
        return result;
    }
    /**
     * Random xorshift, mcg multiply, fixed xorshift
     *
     * The most statistically powerful generator, but all those steps
     * make it slower than some of the others.  We give it the rottenest jobs.
     *
     * Because it's usually used in contexts where the state type and the
     * result type are the same, it is a permutation and is thus invertable.
     * We thus provide a function to invert it.  This function is used to
     * for the "inside out" generator used by the extended generator.
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::rxs_m_xs(WorkI internal) noexcept -> OutI {
        constexpr uint8_t Opbits = [=]() noexcept {
            if (OutBits >= 128)
                return 6;
            if (OutBits >= 64)
                return 5;
            if (OutBits >= 32)
                return 4;
            if (OutBits >= 16)
                return 3;
            return 2;
        }();
        constexpr uint8_t Shift = WorkBits - OutBits;
        constexpr uint8_t Mask = (1 << Opbits) - 1;
        uint8_t rshift = Opbits ? uint8_t(internal >> (WorkBits - Opbits)) & Mask : 0;
        internal ^= internal >> (Opbits + rshift);
        internal *= multiplierMCG();
        OutI result = internal >> Shift;
        result ^= result >> ((2U * OutBits + 2U) / 3U);
        return result;
    }
    /*
     * XSL -- fixed xorshift (to low bits)
     *
     * You shouldn't use this at 64-bits or less.
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::xsl(WorkI internal) noexcept -> OutI {
        constexpr uint8_t TopSpare = SpareBits;
        constexpr uint8_t BottomSpare = SpareBits - TopSpare;
        constexpr uint8_t xshift = (TopSpare + OutBits) / 2;

        internal ^= internal >> xshift;
        OutI result = internal >> BottomSpare;
        return result;
    }
    /*
     * XSL RR -- fixed xorshift (to low bits), random rotate
     *
     * Useful for 128-bit types that are split across two CPU registers.
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::xsl_rr(WorkI internal) noexcept -> OutI {
        constexpr uint8_t WantedOpbits = [=]() noexcept {
            if (OutBits >= 128)
                return 7;
            if (OutBits >= 64)
                return 6;
            if (OutBits >= 32)
                return 5;
            if (OutBits >= 16)
                return 4;
            return 3;
        }();
        constexpr uint8_t Opbits = SpareBits >= WantedOpbits ? WantedOpbits : SpareBits;
        constexpr uint8_t amplifier = WantedOpbits - Opbits;
        constexpr uint8_t Mask = (1 << Opbits) - 1;
        constexpr uint8_t TopSpare = SpareBits;
        constexpr uint8_t BottomSpare = SpareBits - TopSpare;
        constexpr uint8_t xshift = (TopSpare + OutBits) / 2;

        uint8_t rot = Opbits ? uint8_t(internal >> (WorkBits - Opbits)) & Mask : 0;
        uint8_t amprot = (rot << amplifier) & Mask;
        internal ^= internal >> xshift;
        OutI result = OutI(internal >> BottomSpare);
        result = std::rotr(result, amprot);
        return result;
    }
    /*
     * XSL RR RR -- fixed xorshift (to low bits), random rotate (both parts)
     *
     * Useful for 128-bit types that are split across two CPU registers.
     * If you really want an invertable 128-bit RNG, I guess this is the one.
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::xsl_rr_rr(WorkI internal) noexcept -> OutI {
        constexpr uint8_t WantedOpbits = [=]() noexcept {
            if (HalfBits >= 128)
                return 7;
            if (HalfBits >= 64)
                return 6;
            if (HalfBits >= 32)
                return 5;
            if (HalfBits >= 16)
                return 4;
            return 3;
        }();

        constexpr uint8_t Opbits = SpareBits >= WantedOpbits ? WantedOpbits : SpareBits;
        constexpr uint8_t amplifier = WantedOpbits - Opbits;
        constexpr uint8_t Mask = (1 << Opbits) - 1;
        constexpr uint8_t TopSpare = SpareBits;
        constexpr uint8_t xshift = (TopSpare + HalfBits) / 2;

        uint8_t rot = Opbits ? uint8_t(internal >> (WorkBits - Opbits)) & Mask : 0;
        uint8_t amprot = (rot << amplifier) & Mask;
        internal ^= internal >> xshift;
        HalfI lowbits = internal;
        lowbits = std::rotr(lowbits, amprot);
        HalfI highbits = internal >> TopSpare;
        uint8_t rot2 = lowbits & Mask;
        uint8_t amprot2 = (rot2 << amplifier) & Mask;
        highbits = std::rotr(highbits, amprot2);
        return (WorkI(highbits) << TopSpare) ^ WorkI(lowbits);
    }
    /*
     * DXSM -- double xorshift multiply
     *
     * This is a new, more powerful output permutation (added in 2019).  It's
     * a more comprehensive scrambling than RXS_M, but runs faster on 128-bit
     * types.  Although primarily intended for use at large sizes, also works
     * at smaller sizes as well.
     *
     * This permutation is similar to xorshift multiply hash functions, except
     * that one of the multipliers is the LCG multiplier (to avoid needing to
     * have a second constant) and the other is based on the low-order bits.
     * This latter aspect means that the scrambling applied to the high bits
     * depends on the low bits, and makes it (to my eye) impractical to back
     * out the permutation without having the low-order bits.
     */
    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::dxsm(WorkI internal) noexcept -> OutI {
        OutI hi = OutI(internal >> SpareBits);
        OutI lo = OutI(internal);

        lo |= 1;
        hi ^= hi >> (OutBits / 2);
        hi *= OutI(multiplier());
        hi ^= hi >> (3 * (OutBits / 4));
        hi *= lo;
        return hi;
    }

    template<class OutI, class WorkI, OutputFnType Type>
    auto OutputFn<OutI, WorkI, Type>::rxs_m_xs_inv(WorkI internal) noexcept -> WorkI {
        constexpr uint8_t Opbits = [=]() noexcept {
            if (WorkBits >= 128)
                return 6;
            if (WorkBits >= 64)
                return 5;
            if (WorkBits >= 32)
                return 4;
            if (WorkBits >= 16)
                return 3;
            return 2;
        }();
        constexpr uint8_t Mask = (1 << Opbits) - 1;

        internal = unxorshift(internal, WorkBits, (2U * WorkBits + 2U) / 3U);
        internal *= unmultiplierMCG();

        uint8_t rshift = Opbits ? (internal >> (WorkBits - Opbits)) & Mask : 0;
        internal = unxorshift(internal, WorkBits, Opbits + rshift);
        return internal;
    }
}
