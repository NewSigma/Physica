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
/**
 * PCG generator is adapted from [1]
 *
 * Reference:
 * [1] Melissa O'Neill, pcg-cpp (2022); https://github.com/imneme/pcg-cpp
 * [2] PCG generator; http://www.pcg-random.org
 */
#include "PCGImpl/Engine.h"
#include "PCGImpl/Extend.h"
#include "PCGImpl/Stream.h"

namespace Physica::Internal {
    template<typename xtype, typename itype, template<typename XT, typename IT> class output_mixin, bool output_previous = (sizeof(itype) <= 8)>
    using oneseq_base = engine<xtype, itype, output_mixin<xtype, itype>, output_previous, oneseq_stream<itype>>;

    template<typename xtype, typename itype, template<typename XT, typename IT> class output_mixin, bool output_previous = (sizeof(itype) <= 8)>
    using unique_base = engine<xtype, itype, output_mixin<xtype, itype>, output_previous, unique_stream<itype>>;

    template<typename xtype, typename itype, template<typename XT, typename IT> class output_mixin, bool output_previous = (sizeof(itype) <= 8)>
    using setseq_base = engine<xtype, itype, output_mixin<xtype, itype>, output_previous, specific_stream<itype>>;

    template<typename xtype, typename itype, template<typename XT, typename IT> class output_mixin, bool output_previous = (sizeof(itype) <= 8)>
    using mcg_base = engine<xtype, itype, output_mixin<xtype, itype>, output_previous, no_stream<itype>>;

    namespace pcg_engines {
        /* Predefined types for XSH RS */

        using oneseq_xsh_rs_16_8 = oneseq_base<uint8_t, uint16_t, xsh_rs_mixin>;
        using oneseq_xsh_rs_32_16 = oneseq_base<uint16_t, uint32_t, xsh_rs_mixin>;
        using oneseq_xsh_rs_64_32 = oneseq_base<uint32_t, uint64_t, xsh_rs_mixin>;

        using unique_xsh_rs_16_8 = unique_base<uint8_t, uint16_t, xsh_rs_mixin>;
        using unique_xsh_rs_32_16 = unique_base<uint16_t, uint32_t, xsh_rs_mixin>;
        using unique_xsh_rs_64_32 = unique_base<uint32_t, uint64_t, xsh_rs_mixin>;

        using setseq_xsh_rs_16_8 = setseq_base<uint8_t, uint16_t, xsh_rs_mixin>;
        using setseq_xsh_rs_32_16 = setseq_base<uint16_t, uint32_t, xsh_rs_mixin>;
        using setseq_xsh_rs_64_32 = setseq_base<uint32_t, uint64_t, xsh_rs_mixin>;

        using mcg_xsh_rs_16_8 = mcg_base<uint8_t, uint16_t, xsh_rs_mixin>;
        using mcg_xsh_rs_32_16 = mcg_base<uint16_t, uint32_t, xsh_rs_mixin>;
        using mcg_xsh_rs_64_32 = mcg_base<uint32_t, uint64_t, xsh_rs_mixin>;

        /* Predefined types for XSH RR */

        using oneseq_xsh_rr_16_8 = oneseq_base<uint8_t, uint16_t, xsh_rr_mixin>;
        using oneseq_xsh_rr_32_16 = oneseq_base<uint16_t, uint32_t, xsh_rr_mixin>;
        using oneseq_xsh_rr_64_32 = oneseq_base<uint32_t, uint64_t, xsh_rr_mixin>;

        using unique_xsh_rr_16_8 = unique_base<uint8_t, uint16_t, xsh_rr_mixin>;
        using unique_xsh_rr_32_16 = unique_base<uint16_t, uint32_t, xsh_rr_mixin>;
        using unique_xsh_rr_64_32 = unique_base<uint32_t, uint64_t, xsh_rr_mixin>;

        using setseq_xsh_rr_16_8 = setseq_base<uint8_t, uint16_t, xsh_rr_mixin>;
        using setseq_xsh_rr_32_16 = setseq_base<uint16_t, uint32_t, xsh_rr_mixin>;
        using setseq_xsh_rr_64_32 = setseq_base<uint32_t, uint64_t, xsh_rr_mixin>;

        using mcg_xsh_rr_16_8 = mcg_base<uint8_t, uint16_t, xsh_rr_mixin>;
        using mcg_xsh_rr_32_16 = mcg_base<uint16_t, uint32_t, xsh_rr_mixin>;
        using mcg_xsh_rr_64_32 = mcg_base<uint32_t, uint64_t, xsh_rr_mixin>;

        /* Predefined types for RXS M XS */

        using oneseq_rxs_m_xs_8_8 = oneseq_base<uint8_t, uint8_t, rxs_m_xs_mixin>;
        using oneseq_rxs_m_xs_16_16 = oneseq_base<uint16_t, uint16_t, rxs_m_xs_mixin>;
        using oneseq_rxs_m_xs_32_32 = oneseq_base<uint32_t, uint32_t, rxs_m_xs_mixin>;
        using oneseq_rxs_m_xs_64_64 = oneseq_base<uint64_t, uint64_t, rxs_m_xs_mixin>;

        using unique_rxs_m_xs_8_8 = unique_base<uint8_t, uint8_t, rxs_m_xs_mixin>;
        using unique_rxs_m_xs_16_16 = unique_base<uint16_t, uint16_t, rxs_m_xs_mixin>;
        using unique_rxs_m_xs_32_32 = unique_base<uint32_t, uint32_t, rxs_m_xs_mixin>;
        using unique_rxs_m_xs_64_64 = unique_base<uint64_t, uint64_t, rxs_m_xs_mixin>;

        using setseq_rxs_m_xs_8_8 = setseq_base<uint8_t, uint8_t, rxs_m_xs_mixin>;
        using setseq_rxs_m_xs_16_16 = setseq_base<uint16_t, uint16_t, rxs_m_xs_mixin>;
        using setseq_rxs_m_xs_32_32 = setseq_base<uint32_t, uint32_t, rxs_m_xs_mixin>;
        using setseq_rxs_m_xs_64_64 = setseq_base<uint64_t, uint64_t, rxs_m_xs_mixin>;

        // MCG versions don't make sense here, so aren't defined.

        /* Predefined types for RXS M */

        using oneseq_rxs_m_16_8 = oneseq_base<uint8_t, uint16_t, rxs_m_mixin>;
        using oneseq_rxs_m_32_16 = oneseq_base<uint16_t, uint32_t, rxs_m_mixin>;
        using oneseq_rxs_m_64_32 = oneseq_base<uint32_t, uint64_t, rxs_m_mixin>;

        using unique_rxs_m_16_8 = unique_base<uint8_t, uint16_t, rxs_m_mixin>;
        using unique_rxs_m_32_16 = unique_base<uint16_t, uint32_t, rxs_m_mixin>;
        using unique_rxs_m_64_32 = unique_base<uint32_t, uint64_t, rxs_m_mixin>;

        using setseq_rxs_m_16_8 = setseq_base<uint8_t, uint16_t, rxs_m_mixin>;
        using setseq_rxs_m_32_16 = setseq_base<uint16_t, uint32_t, rxs_m_mixin>;
        using setseq_rxs_m_64_32 = setseq_base<uint32_t, uint64_t, rxs_m_mixin>;

        using mcg_rxs_m_16_8 = mcg_base<uint8_t, uint16_t, rxs_m_mixin>;
        using mcg_rxs_m_32_16 = mcg_base<uint16_t, uint32_t, rxs_m_mixin>;
        using mcg_rxs_m_64_32 = mcg_base<uint32_t, uint64_t, rxs_m_mixin>;

        /* Predefined types for DXSM */

        using oneseq_dxsm_16_8 = oneseq_base<uint8_t, uint16_t, dxsm_mixin>;
        using oneseq_dxsm_32_16 = oneseq_base<uint16_t, uint32_t, dxsm_mixin>;
        using oneseq_dxsm_64_32 = oneseq_base<uint32_t, uint64_t, dxsm_mixin>;

        using unique_dxsm_16_8 = unique_base<uint8_t, uint16_t, dxsm_mixin>;
        using unique_dxsm_32_16 = unique_base<uint16_t, uint32_t, dxsm_mixin>;
        using unique_dxsm_64_32 = unique_base<uint32_t, uint64_t, dxsm_mixin>;

        using setseq_dxsm_16_8 = setseq_base<uint8_t, uint16_t, dxsm_mixin>;
        using setseq_dxsm_32_16 = setseq_base<uint16_t, uint32_t, dxsm_mixin>;
        using setseq_dxsm_64_32 = setseq_base<uint32_t, uint64_t, dxsm_mixin>;

        using mcg_dxsm_16_8 = mcg_base<uint8_t, uint16_t, dxsm_mixin>;
        using mcg_dxsm_32_16 = mcg_base<uint16_t, uint32_t, dxsm_mixin>;
        using mcg_dxsm_64_32 = mcg_base<uint32_t, uint64_t, dxsm_mixin>;

        /* Predefined types for XSL RR (only defined for "large" types) */

        using oneseq_xsl_rr_64_32 = oneseq_base<uint32_t, uint64_t, xsl_rr_mixin>;

        using unique_xsl_rr_64_32 = unique_base<uint32_t, uint64_t, xsl_rr_mixin>;

        using setseq_xsl_rr_64_32 = setseq_base<uint32_t, uint64_t, xsl_rr_mixin>;

        using mcg_xsl_rr_64_32 = mcg_base<uint32_t, uint64_t, xsl_rr_mixin>;

        /* Predefined types for XSL RR RR (only defined for "large" types) */

        using oneseq_xsl_rr_rr_64_64 = oneseq_base<uint64_t, uint64_t, xsl_rr_rr_mixin>;

        using unique_xsl_rr_rr_64_64 = unique_base<uint64_t, uint64_t, xsl_rr_rr_mixin>;

        using setseq_xsl_rr_rr_64_64 = setseq_base<uint64_t, uint64_t, xsl_rr_rr_mixin>;

        // MCG versions don't make sense here, so aren't defined.

        /* Extended generators */

        template<uint8_t table_pow2, uint8_t advance_pow2, typename BaseRNG, bool KDD = true>
        using ext_std8 = extended<table_pow2, advance_pow2, BaseRNG, oneseq_rxs_m_xs_8_8, KDD>;

        template<uint8_t table_pow2, uint8_t advance_pow2, typename BaseRNG, bool KDD = true>
        using ext_std16 = extended<table_pow2, advance_pow2, BaseRNG, oneseq_rxs_m_xs_16_16, KDD>;

        template<uint8_t table_pow2, uint8_t advance_pow2, typename BaseRNG, bool KDD = true>
        using ext_std32 = extended<table_pow2, advance_pow2, BaseRNG, oneseq_rxs_m_xs_32_32, KDD>;

        template<uint8_t table_pow2, uint8_t advance_pow2, typename BaseRNG, bool KDD = true>
        using ext_std64 = extended<table_pow2, advance_pow2, BaseRNG, oneseq_rxs_m_xs_64_64, KDD>;

        template<uint8_t table_pow2, uint8_t advance_pow2, bool KDD = true>
        using ext_oneseq_rxs_m_xs_32_32 =
                ext_std32<table_pow2, advance_pow2, oneseq_rxs_m_xs_32_32, KDD>;

        template<uint8_t table_pow2, uint8_t advance_pow2, bool KDD = true>
        using ext_mcg_xsh_rs_64_32 =
                ext_std32<table_pow2, advance_pow2, mcg_xsh_rs_64_32, KDD>;

        template<uint8_t table_pow2, uint8_t advance_pow2, bool KDD = true>
        using ext_oneseq_xsh_rs_64_32 =
                ext_std32<table_pow2, advance_pow2, oneseq_xsh_rs_64_32, KDD>;

        template<uint8_t table_pow2, uint8_t advance_pow2, bool KDD = true>
        using ext_setseq_xsh_rr_64_32 = ext_std32<table_pow2, advance_pow2, setseq_xsh_rr_64_32, KDD>;
    }

    using pcg32 = pcg_engines::setseq_xsh_rr_64_32;
    using pcg32_oneseq = pcg_engines::oneseq_xsh_rr_64_32;
    using pcg32_unique = pcg_engines::unique_xsh_rr_64_32;
    using pcg32_fast = pcg_engines::mcg_xsh_rs_64_32;

    using pcg8_once_insecure = pcg_engines::setseq_rxs_m_xs_8_8;
    using pcg16_once_insecure = pcg_engines::setseq_rxs_m_xs_16_16;
    using pcg32_once_insecure = pcg_engines::setseq_rxs_m_xs_32_32;
    using pcg64_once_insecure = pcg_engines::setseq_rxs_m_xs_64_64;

    using pcg8_oneseq_once_insecure = pcg_engines::oneseq_rxs_m_xs_8_8;
    using pcg16_oneseq_once_insecure = pcg_engines::oneseq_rxs_m_xs_16_16;
    using pcg32_oneseq_once_insecure = pcg_engines::oneseq_rxs_m_xs_32_32;
    using pcg64_oneseq_once_insecure = pcg_engines::oneseq_rxs_m_xs_64_64;

    // These two extended RNGs provide two-dimensionally equidistributed
    // 32-bit generators.  pcg32_k2_fast occupies the same space as pcg64,
    // and can be called twice to generate 64 bits, but does not required
    // 128-bit math; on 32-bit systems, it's faster than pcg64 as well.

    using pcg32_k2 = pcg_engines::ext_setseq_xsh_rr_64_32<1, 16, true>;
    using pcg32_k2_fast = pcg_engines::ext_oneseq_xsh_rs_64_32<1, 32, true>;

    // These eight extended RNGs have about as much state as arc4random
    //
    //  - the k variants are k-dimensionally equidistributed
    //  - the c variants offer are intended to be harder to predict
    //
    // (neither is intended for use in cryptographic applications)

    using pcg32_k64 = pcg_engines::ext_setseq_xsh_rr_64_32<6, 16, true>;
    using pcg32_k64_oneseq = pcg_engines::ext_mcg_xsh_rs_64_32<6, 32, true>;
    using pcg32_k64_fast = pcg_engines::ext_oneseq_xsh_rs_64_32<6, 32, true>;

    using pcg32_c64 = pcg_engines::ext_setseq_xsh_rr_64_32<6, 16, false>;
    using pcg32_c64_oneseq = pcg_engines::ext_oneseq_xsh_rs_64_32<6, 32, false>;
    using pcg32_c64_fast = pcg_engines::ext_mcg_xsh_rs_64_32<6, 32, false>;

    // These eight extended RNGs have more state than the Mersenne twister
    //
    //  - the k variants are k-dimensionally equidistributed
    //  - the c variants offer are intended to be harder to predict
    //
    // (neither is intended for use in cryptographic applications)

    using pcg32_k1024 = pcg_engines::ext_setseq_xsh_rr_64_32<10, 16, true>;
    using pcg32_k1024_fast = pcg_engines::ext_oneseq_xsh_rs_64_32<10, 32, true>;

    using pcg32_c1024 = pcg_engines::ext_setseq_xsh_rr_64_32<10, 16, false>;
    using pcg32_c1024_fast = pcg_engines::ext_oneseq_xsh_rs_64_32<10, 32, false>;

    // These generators have an insanely huge period (2^524352), and is suitable
    // for silly party tricks, such as dumping out 64 KB ZIP files at an arbitrary
    // point in the future.   [Actually, over the full period of the generator, it
    // will produce every 64 KB ZIP file 2^64 times!]

    using pcg32_k16384 = pcg_engines::ext_setseq_xsh_rr_64_32<14, 16, true>;
    using pcg32_k16384_fast = pcg_engines::ext_oneseq_xsh_rs_64_32<14, 32, true>;
}
