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

namespace Physica::Internal {
    template<class T>
    consteval static T multiplier() {
        if constexpr (std::same_as<uint8_t, T>)
            return 141U;
        if constexpr (std::same_as<uint16_t, T>)
            return 12829U;
        if constexpr (std::same_as<uint32_t, T>)
            return 747796405U;
        if constexpr (std::same_as<uint64_t, T>)
            return 6364136223846793005ULL;
        else {
            static_assert(std::same_as<__uint128_t, T>);
            return (__uint128_t(2549297995355413924ULL) << 64U) + 4865540595714422341ULL;
        }
    }

    template<class T>
    consteval static T increment() {
        if constexpr (std::same_as<uint8_t, T>)
            return 77U;
        if constexpr (std::same_as<uint16_t, T>)
            return 47989U;
        if constexpr (std::same_as<uint32_t, T>)
            return 2891336453U;
        if constexpr (std::same_as<uint64_t, T>)
            return 1442695040888963407ULL;
        else {
            static_assert(std::same_as<__uint128_t, T>);
            return (__uint128_t(6364136223846793005ULL) << 64U) + 1442695040888963407ULL;
        }
    }

    template<class T>
    consteval static T multiplierMCG() {
        if constexpr (std::same_as<uint8_t, T>)
            return 217U;
        if constexpr (std::same_as<uint16_t, T>)
            return 62169U;
        if constexpr (std::same_as<uint32_t, T>)
            return 277803737U;
        if constexpr (std::same_as<uint64_t, T>)
            return 12605985483714917081ULL;
        else {
            static_assert(std::same_as<__uint128_t, T>);
            return (__uint128_t(17766728186571221404ULL) << 64U) + 12605985483714917081ULL;
        }
    }

    template<class T>
    consteval static T unmultiplierMCG() {
        if constexpr (std::same_as<uint8_t, T>)
            return 105U;
        if constexpr (std::same_as<uint16_t, T>)
            return 28009U;
        if constexpr (std::same_as<uint32_t, T>)
            return 2897767785U;
        if constexpr (std::same_as<uint64_t, T>)
            return 15009553638781119849ULL;
        else {
            static_assert(std::same_as<__uint128_t, T>);
            return (__uint128_t(14422606686972528997ULL) << 64U) + 15009553638781119849ULL;
        }
    }

    /*
     * OUTPUT FUNCTIONS.
     *
     * These are the core of the PCG generation scheme.  They specify how to
     * turn the base LCG's internal state into the output value of the final
     * generator.
     *
     * They're implemented as mixin classes.
     *
     * All of the classes have code that is written to allow it to be applied
     * at *arbitrary* bit sizes, although in practice they'll only be used at
     * standard sizes supported by C++.
     */

    /*
     * XSH RS -- high xorshift, followed by a random shift
     *
     * Fast.  A good performer.
     */
    template<typename xtype, typename itype>
    struct xsh_rs_mixin {
        static xtype output(itype internal) {
            constexpr uint8_t bits = sizeof(itype) * 8;
            constexpr uint8_t xtypebits = sizeof(xtype) * 8;
            constexpr uint8_t sparebits = bits - xtypebits;
            constexpr uint8_t opbits = [=]() {
                if (sparebits >= 64 + 5)
                    return 5;
                if (sparebits >= 32 + 4)
                    return 4;
                if (sparebits >= 16 + 3)
                    return 3;
                if (sparebits >= 4 + 2)
                    return 2;
                return sparebits >= (1 + 1) ? 1 : 0;
            }();
            constexpr uint8_t mask = (1 << opbits) - 1;
            constexpr uint8_t maxrandshift = mask;
            constexpr uint8_t topspare = opbits;
            constexpr uint8_t bottomspare = sparebits - topspare;
            constexpr uint8_t xshift = topspare + (xtypebits + maxrandshift) / 2;
            uint8_t rshift =
                    opbits ? uint8_t(internal >> (bits - opbits)) & mask : 0;
            internal ^= internal >> xshift;
            xtype result = xtype(internal >> (bottomspare - maxrandshift + rshift));
            return result;
        }
    };

    /*
     * XSH RR -- high xorshift, followed by a random rotate
     *
     * Fast.  A good performer.  Slightly better statistically than XSH RS.
     */
    template<typename xtype, typename itype>
    struct xsh_rr_mixin {
        static xtype output(itype internal) {
            constexpr uint8_t bits = uint8_t(sizeof(itype) * 8);
            constexpr uint8_t xtypebits = uint8_t(sizeof(xtype) * 8);
            constexpr uint8_t sparebits = bits - xtypebits;
            constexpr uint8_t wantedopbits = [=]() {
                if (xtypebits >= 128)
                    return 7;
                if (xtypebits >= 64)
                    return 6;
                if (xtypebits >= 32)
                    return 5;
                if (xtypebits >= 16)
                    return 4;
                return 3;
            }();
            constexpr uint8_t opbits = sparebits >= wantedopbits ? wantedopbits : sparebits;
            constexpr uint8_t amplifier = wantedopbits - opbits;
            constexpr uint8_t mask = (1 << opbits) - 1;
            constexpr uint8_t topspare = opbits;
            constexpr uint8_t bottomspare = sparebits - topspare;
            constexpr uint8_t xshift = (topspare + xtypebits) / 2;
            uint8_t rot = opbits ? uint8_t(internal >> (bits - opbits)) & mask
                                 : 0;
            uint8_t amprot = (rot << amplifier) & mask;
            internal ^= internal >> xshift;
            xtype result = xtype(internal >> bottomspare);
            result = rotr(result, amprot);
            return result;
        }
    };

    /*
     * RXS -- random xorshift
     */
    template<typename xtype, typename itype>
    struct rxs_mixin {
        static xtype output_rxs(itype internal) {
            constexpr uint8_t bits = sizeof(itype) * 8;
            constexpr uint8_t xtypebits = sizeof(xtype) * 8;
            constexpr uint8_t shift = bits - xtypebits;
            constexpr uint8_t extrashift = (xtypebits - shift) / 2;
            uint8_t rshift = [=]() {
                if (shift > (64 + 8))
                    return (internal >> (bits - 6)) & 63;
                if (shift > (32 + 4))
                    return (internal >> (bits - 5)) & 31;
                if (shift > (16 + 2))
                    return (internal >> (bits - 4)) & 15;
                if (shift > (8 + 1))
                    return (internal >> (bits - 3)) & 7;
                if (shift > (4 + 1))
                    return (internal >> (bits - 2)) & 3;
                return (shift > (2 + 1)) ? (internal >> (bits - 1)) & 1 : 0;
            }();

            internal ^= internal >> (shift + extrashift - rshift);
            xtype result = internal >> rshift;
            return result;
        }
    };

    /*
     * RXS M XS -- random xorshift, mcg multiply, fixed xorshift
     *
     * The most statistically powerful generator, but all those steps
     * make it slower than some of the others.  We give it the rottenest jobs.
     *
     * Because it's usually used in contexts where the state type and the
     * result type are the same, it is a permutation and is thus invertable.
     * We thus provide a function to invert it.  This function is used to
     * for the "inside out" generator used by the extended generator.
     */
    template<typename xtype, typename itype>
    struct rxs_m_xs_mixin {
        static xtype output(itype internal) {
            constexpr uint8_t xtypebits = sizeof(xtype) * 8;
            constexpr uint8_t bits = sizeof(itype) * 8;
            constexpr uint8_t opbits = [=]() {
                if (xtypebits >= 128)
                    return 6;
                if (xtypebits >= 64)
                    return 5;
                if (xtypebits >= 32)
                    return 4;
                if (xtypebits >= 16)
                    return 3;
                return 2;
            }();
            constexpr uint8_t shift = bits - xtypebits;
            constexpr uint8_t mask = (1 << opbits) - 1;
            uint8_t rshift = opbits ? uint8_t(internal >> (bits - opbits)) & mask : 0;
            internal ^= internal >> (opbits + rshift);
            internal *= multiplierMCG<itype>();
            xtype result = internal >> shift;
            result ^= result >> ((2U * xtypebits + 2U) / 3U);
            return result;
        }

        static itype unoutput(itype internal) {
            constexpr uint8_t bits = sizeof(itype) * 8;
            constexpr uint8_t opbits = [=]() {
                if (bits >= 128)
                    return 6;
                if (bits >= 64)
                    return 5;
                if (bits >= 32)
                    return 4;
                if (bits >= 16)
                    return 3;
                return 2;
            }();
            constexpr uint8_t mask = (1 << opbits) - 1;

            internal = unxorshift(internal, bits, (2U * bits + 2U) / 3U);
            internal *= unmultiplierMCG<itype>();

            uint8_t rshift = opbits ? (internal >> (bits - opbits)) & mask : 0;
            internal = unxorshift(internal, bits, opbits + rshift);
            return internal;
        }
    };

    /*
     * RXS M -- random xorshift, mcg multiply
     */
    template<typename xtype, typename itype>
    struct rxs_m_mixin {
        static xtype output(itype internal) {
            constexpr uint8_t xtypebits = sizeof(xtype) * 8;
            constexpr uint8_t bits = sizeof(itype) * 8;
            constexpr uint8_t opbits = [=]() {
                if (xtypebits >= 128)
                    return 6;
                if (xtypebits >= 64)
                    return 5;
                if (xtypebits >= 32)
                    return 4;
                if (xtypebits >= 16)
                    return 3;
                return 2;
            }();
            constexpr uint8_t shift = bits - xtypebits;
            constexpr uint8_t mask = (1 << opbits) - 1;
            uint8_t rshift = opbits ? (internal >> (bits - opbits)) & mask : 0;
            internal ^= internal >> (opbits + rshift);
            internal *= multiplierMCG<itype>();
            xtype result = internal >> shift;
            return result;
        }
    };

    /*
     * DXSM -- double xorshift multiply
     *
     * This is a new, more powerful output permutation (added in 2019).  It's
     * a more comprehensive scrambling than RXS M, but runs faster on 128-bit
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
    template<typename xtype, typename itype>
    struct dxsm_mixin {
        xtype output(itype internal) {
            constexpr uint8_t xtypebits = sizeof(xtype) * 8;
            constexpr uint8_t itypebits = sizeof(itype) * 8;
            static_assert(xtypebits <= itypebits / 2,
                          "Output type must be half the size of the state type.");

            xtype hi = xtype(internal >> (itypebits - xtypebits));
            xtype lo = xtype(internal);

            lo |= 1;
            hi ^= hi >> (xtypebits / 2);
            hi *= xtype(multiplier<itype>());
            hi ^= hi >> (3 * (xtypebits / 4));
            hi *= lo;
            return hi;
        }
    };

    /*
     * XSL RR -- fixed xorshift (to low bits), random rotate
     *
     * Useful for 128-bit types that are split across two CPU registers.
     */
    template<typename xtype, typename itype>
    struct xsl_rr_mixin {
        static xtype output(itype internal) {
            constexpr uint8_t xtypebits = sizeof(xtype) * 8;
            constexpr uint8_t bits = sizeof(itype) * 8;
            constexpr uint8_t sparebits = bits - xtypebits;
            constexpr uint8_t wantedopbits = [=]() {
                if (xtypebits >= 128)
                    return 7;
                if (xtypebits >= 64)
                    return 6;
                if (xtypebits >= 32)
                    return 5;
                if (xtypebits >= 16)
                    return 4;
                return 3;
            }();
            constexpr uint8_t opbits = sparebits >= wantedopbits ? wantedopbits : sparebits;
            constexpr uint8_t amplifier = wantedopbits - opbits;
            constexpr uint8_t mask = (1 << opbits) - 1;
            constexpr uint8_t topspare = sparebits;
            constexpr uint8_t bottomspare = sparebits - topspare;
            constexpr uint8_t xshift = (topspare + xtypebits) / 2;

            uint8_t rot =
                    opbits ? uint8_t(internal >> (bits - opbits)) & mask : 0;
            uint8_t amprot = (rot << amplifier) & mask;
            internal ^= internal >> xshift;
            xtype result = xtype(internal >> bottomspare);
            result = rotr(result, amprot);
            return result;
        }
    };

    /*
     * XSL RR RR -- fixed xorshift (to low bits), random rotate (both parts)
     *
     * Useful for 128-bit types that are split across two CPU registers.
     * If you really want an invertable 128-bit RNG, I guess this is the one.
     */
    template<typename xtype, typename itype>
    struct xsl_rr_rr_mixin {
        template<typename T> struct halfsize_trait;
        template<> struct halfsize_trait<uint64_t> {
            using type = uint32_t;
        };
        template<> struct halfsize_trait<uint32_t> {
            using type = uint16_t;
        };
        template<> struct halfsize_trait<uint16_t> {
            using type = uint8_t;
        };
        using htype = typename halfsize_trait<itype>::type;

        static itype output(itype internal) {
            constexpr uint8_t htypebits = sizeof(htype) * 8;
            constexpr uint8_t bits = sizeof(itype) * 8;
            constexpr uint8_t sparebits = bits - htypebits;
            constexpr uint8_t wantedopbits = [=]() {
                if (htypebits >= 128)
                    return 7;
                if (htypebits >= 64)
                    return 6;
                if (htypebits >= 32)
                    return 5;
                if (htypebits >= 16)
                    return 4;
                return 3;
            }();

            constexpr uint8_t opbits = sparebits >= wantedopbits ? wantedopbits : sparebits;
            constexpr uint8_t amplifier = wantedopbits - opbits;
            constexpr uint8_t mask = (1 << opbits) - 1;
            constexpr uint8_t topspare = sparebits;
            constexpr uint8_t xshift = (topspare + htypebits) / 2;

            uint8_t rot =
                    opbits ? uint8_t(internal >> (bits - opbits)) & mask : 0;
            uint8_t amprot = (rot << amplifier) & mask;
            internal ^= internal >> xshift;
            htype lowbits = htype(internal);
            lowbits = rotr(lowbits, amprot);
            htype highbits = htype(internal >> topspare);
            uint8_t rot2 = lowbits & mask;
            uint8_t amprot2 = (rot2 << amplifier) & mask;
            highbits = rotr(highbits, amprot2);
            return (itype(highbits) << topspare) ^ itype(lowbits);
        }
    };

    /*
     * XSH -- fixed xorshift (to high bits)
     *
     * You shouldn't use this at 64-bits or less.
     */
    template<typename xtype, typename itype>
    struct xsh_mixin {
        static xtype output(itype internal) {
            constexpr uint8_t xtypebits = sizeof(xtype) * 8;
            constexpr uint8_t bits = sizeof(itype) * 8;
            constexpr uint8_t sparebits = bits - xtypebits;
            constexpr uint8_t topspare = 0;
            constexpr uint8_t bottomspare = sparebits - topspare;
            constexpr uint8_t xshift = (topspare + xtypebits) / 2;

            internal ^= internal >> xshift;
            xtype result = internal >> bottomspare;
            return result;
        }
    };

    /*
     * XSL -- fixed xorshift (to low bits)
     *
     * You shouldn't use this at 64-bits or less.
     */
    template<typename xtype, typename itype>
    struct xsl_mixin {
        xtype output(itype internal) {
            constexpr uint8_t xtypebits = sizeof(xtype) * 8;
            constexpr uint8_t bits = sizeof(itype) * 8;
            constexpr uint8_t sparebits = bits - xtypebits;
            constexpr uint8_t topspare = sparebits;
            constexpr uint8_t bottomspare = sparebits - topspare;
            constexpr uint8_t xshift = (topspare + xtypebits) / 2;

            internal ^= internal >> xshift;
            xtype result = internal >> bottomspare;
            return result;
        }
    };
}
