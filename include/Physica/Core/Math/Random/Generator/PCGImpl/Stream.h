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

#include <cstddef>
#include "Physica/Core/Exception/NoImplException.h"
#include "OutputFn.h"

namespace Physica::Internal {
    /*
     * Each PCG generator is available in four variants, based on how it applies
     * the additive constant for its underlying LCG; the variations are:
     *
     *     single stream   - all instances use the same fixed constant, thus
     *                       the RNG always somewhere in same sequence
     *     mcg             - adds zero, resulting in a single stream and reduced
     *                       period
     *     specific stream - the constant can be changed at any time, selecting
     *                       a different random sequence
     *     unique stream   - the constant is based on the memory address of the
     *                       object, thus every RNG has its own unique sequence
     *
     * This variation is provided though mixin classes which define a function
     * value called increment() that returns the necessary additive constant.
     */
    template<typename itype>
    class unique_stream {
    protected:
        constexpr static bool is_mcg = false;

        void set_stream(...) {
            Physica::noImpl("Is never called, but is provided for symmetry with specific_stream");
        }
    public:
        using state_type = itype;

        constexpr itype increment() const {
            return itype(reinterpret_cast<uintptr_t>(this) | 1);
        }

        constexpr itype stream() const {
            return increment() >> 1;
        }

        constexpr static bool can_specify_stream = false;

        constexpr static size_t streams_pow2() {
            return (sizeof(itype) < sizeof(size_t) ? sizeof(itype)
                                                   : sizeof(size_t))
                         * 8
                 - 1U;
        }
    protected:
        constexpr unique_stream() = default;
    };

    template<typename itype>
    class no_stream {
    protected:
        constexpr static bool is_mcg = true;

        void set_stream(...) {
            Physica::noImpl("Is never called, but is provided for symmetry with specific_stream");
        }
    public:
        using state_type = itype;

        constexpr static itype increment() {
            return 0;
        }

        constexpr static bool can_specify_stream = false;

        constexpr static size_t streams_pow2() {
            return 0U;
        }
    protected:
        constexpr no_stream() = default;
    };

    template<typename itype>
    class oneseq_stream {
    protected:
        constexpr static bool is_mcg = false;

        void set_stream(...) {
            Physica::noImpl("Is never called, but is provided for symmetry with specific_stream");
        }
    public:
        using state_type = itype;

        constexpr static itype stream() {
            return increment<itype>() >> 1;
        }

        constexpr static bool can_specify_stream = false;

        constexpr static size_t streams_pow2() {
            return 0U;
        }
    protected:
        constexpr oneseq_stream() = default;
    };

    template<typename itype>
    class specific_stream {
    protected:
        constexpr static bool is_mcg = false;

        itype inc_ = Internal::increment<itype>();
    public:
        using state_type = itype;
        using stream_state = itype;

        constexpr itype increment() const {
            return inc_;
        }

        itype stream() {
            return inc_ >> 1;
        }

        void set_stream(itype specific_seq) {
            inc_ = (specific_seq << 1) | 1;
        }

        constexpr static bool can_specify_stream = true;

        constexpr static size_t streams_pow2() {
            return (sizeof(itype) * 8) - 1U;
        }
    protected:
        specific_stream() = default;

        specific_stream(itype specific_seq)
                : inc_(itype(specific_seq << 1) | itype(1U)) {}
    };
}
