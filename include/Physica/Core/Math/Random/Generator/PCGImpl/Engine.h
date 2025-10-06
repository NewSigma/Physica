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
#include <iostream>
#include <type_traits>
#include <utility>
#include "OutputFn.h"

namespace Physica::Internal {
    /*
     * This is where it all comes together.  This function joins together three
     * mixin classes which define
     *    - the LCG additive constant (the stream)
     *    - the LCG multiplier
     *    - the output function
     * in addition, we specify the type of the LCG state, and the result type,
     * and whether to use the pre-advance version of the state for the output
     * (increasing instruction-level parallelism) or the post-advance version
     * (reducing register pressure).
     *
     * Given the high level of parameterization, the code has to use some
     * template-metaprogramming tricks to handle some of the subtle variations
     * involved.
     */
    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    class engine : public stream_mixin, protected output_mixin {
    protected:
        itype state_;

        using stream_mixin::increment;
    public:
        using result_type = xtype;
        using state_type = itype;
    public:
        explicit engine(itype state = itype(0xcafef00dd15ea5e5ULL)) noexcept;
        template<typename sm = stream_mixin>
        engine(itype state, typename sm::stream_state stream_seed) noexcept;
        template<typename SeedSeq>
        explicit engine(SeedSeq& seedSeq) noexcept requires(!stream_mixin::can_specify_stream);
        template<typename SeedSeq>
        explicit engine(SeedSeq& seedSeq) noexcept requires(stream_mixin::can_specify_stream);
        /* Operators */
        result_type operator()() noexcept {
            if (output_previous)
                return this->output(base_generate0());
            return this->output(base_generate());
        }

        result_type operator()(result_type upper_bound) noexcept {
            return bounded_rand(*this, upper_bound);
        }

        template<typename xtype1, typename itype1, typename output_mixin1, bool output_previous1, typename stream_mixin_lhs, typename stream_mixin_rhs>
        friend bool operator==(const engine<xtype1, itype1, output_mixin1, output_previous1, stream_mixin_lhs>&,
                               const engine<xtype1, itype1, output_mixin1, output_previous1, stream_mixin_rhs>&);

        template<typename xtype1, typename itype1, typename output_mixin1, bool output_previous1, typename stream_mixin_lhs, typename stream_mixin_rhs>
        friend itype1 operator-(const engine<xtype1, itype1, output_mixin1, output_previous1, stream_mixin_lhs>&,
                                const engine<xtype1, itype1, output_mixin1, output_previous1, stream_mixin_rhs>&);

        template<typename CharT, typename Traits, typename xtype1, typename itype1, typename output_mixin1, bool output_previous1, typename stream_mixin1>
        friend std::basic_ostream<CharT, Traits>&
        operator<<(std::basic_ostream<CharT, Traits>& out, const engine<xtype1, itype1, output_mixin1, output_previous1, stream_mixin1>&);

        template<typename CharT, typename Traits, typename xtype1, typename itype1, typename output_mixin1, bool output_previous1, typename stream_mixin1>
        friend std::basic_istream<CharT, Traits>&
        operator>>(std::basic_istream<CharT, Traits>& in, engine<xtype1, itype1, output_mixin1, output_previous1, stream_mixin1>& rng);
        /* Operations */
        void advance(itype delta) noexcept;
        void backstep(itype delta) noexcept;
        void discard(itype delta) noexcept;
        bool wrapped() noexcept;

        template<typename... Args>
        void seed(Args&&... args) noexcept;
        /* Static members */
        constexpr static size_t period_pow2() noexcept { return sizeof(state_type) * 8 - 2 * stream_mixin::is_mcg; }
        constexpr static result_type min() noexcept { return std::numeric_limits<result_type>::min(); }
        constexpr static result_type max() noexcept { return std::numeric_limits<result_type>::max(); }
    protected:
        itype bump(itype state) noexcept;
        itype base_generate() noexcept;
        itype base_generate0() noexcept;
        itype distance(itype newstate, itype mask = itype(~itype(0U))) const noexcept;
        /* Static members */
        static itype advance(itype state, itype delta, itype cur_mult, itype cur_plus) noexcept;
        static itype distance(itype cur_state, itype newstate, itype cur_mult, itype cur_plus, itype mask = ~itype(0U)) noexcept;
    };

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    engine<xtype, itype, output_mixin, output_previous, stream_mixin>::engine(itype state) noexcept
            : state_(this->is_mcg ? (state | state_type(3U)) : bump(state + this->increment())) {}

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    template<typename sm>
    engine<xtype, itype, output_mixin, output_previous, stream_mixin>::engine(itype state, typename sm::stream_state stream_seed) noexcept
            : stream_mixin(stream_seed), state_(this->is_mcg ? state | state_type(3U) : bump(state + this->increment())) {}

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    template<typename SeedSeq>
    engine<xtype, itype, output_mixin, output_previous, stream_mixin>::engine(SeedSeq& seedSeq) noexcept requires(!stream_mixin::can_specify_stream)
            : engine(generate_one<itype>(std::forward<SeedSeq>(seedSeq))) {}

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    template<typename SeedSeq>
    engine<xtype, itype, output_mixin, output_previous, stream_mixin>::engine(SeedSeq& seedSeq) noexcept requires(stream_mixin::can_specify_stream)
    {
        std::array<itype, 2> seeddata{};
        seedSeq.generate(seeddata);
        seed(seeddata[1], seeddata[0]);
    }

    template<typename CharT, typename Traits, typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& out, const engine<xtype, itype, output_mixin, output_previous, stream_mixin>& rng) {
        auto orig_flags = out.flags(std::ios_base::dec | std::ios_base::left);
        auto space = out.widen(' ');
        auto orig_fill = out.fill();

        out << rng.multiplier() << space
            << rng.increment() << space
            << rng.state_;

        out.flags(orig_flags);
        out.fill(orig_fill);
        return out;
    }

    template<typename CharT, typename Traits, typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& in, engine<xtype, itype, output_mixin, output_previous, stream_mixin>& rng) {
        auto orig_flags = in.flags(std::ios_base::dec | std::ios_base::skipws);

        itype multiplier;
        itype increment;
        itype state;
        in >> multiplier >> increment >> state;

        bool good = !in.fail() && (multiplier == rng.multiplier()) && (increment == rng.increment());
        if (good) {
            if (rng.can_specify_stream)
                rng.set_stream(increment >> 1);
            rng.state_ = state;
        }
        in.flags(orig_flags);
        return in;
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin_lhs, typename stream_mixin_rhs>
    itype operator-(const engine<xtype, itype, output_mixin, output_previous, stream_mixin_lhs>& lhs,
                    const engine<xtype, itype, output_mixin, output_previous, stream_mixin_rhs>& rhs) {
        static_assert(std::same_as<stream_mixin_lhs, stream_mixin_rhs>, "Incomparable generators");
        if (lhs.increment() == rhs.increment())
            return rhs.distance(lhs.state_);

        constexpr itype ONE = 1U;
        itype lhs_diff = lhs.increment() + (lhs.multiplier() - ONE) * lhs.state_;
        itype rhs_diff = rhs.increment() + (rhs.multiplier() - ONE) * rhs.state_;
        if ((lhs_diff & itype(3U)) != (rhs_diff & itype(3U))) {
            rhs_diff = -rhs_diff;
        }
        return rhs.distance(rhs_diff, lhs_diff, rhs.multiplier(), itype(0U));
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin_lhs, typename stream_mixin_rhs>
    bool operator==(const engine<xtype, itype, output_mixin, output_previous, stream_mixin_lhs>& lhs,
                    const engine<xtype, itype, output_mixin, output_previous, stream_mixin_rhs>& rhs) {
        return (lhs.multiplier() == rhs.multiplier())
            && (lhs.increment() == rhs.increment())
            && (lhs.state_ == rhs.state_);
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin_lhs, typename stream_mixin_rhs>
    inline bool operator!=(const engine<xtype, itype, output_mixin, output_previous, stream_mixin_lhs>& lhs,
                           const engine<xtype, itype, output_mixin, output_previous, stream_mixin_rhs>& rhs) {
        return !operator==(lhs, rhs);
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    void engine<xtype, itype, output_mixin, output_previous, stream_mixin>::advance(itype delta) noexcept {
        state_ = advance(state_, delta, this->multiplier(), this->increment());
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    void engine<xtype, itype, output_mixin, output_previous, stream_mixin>::backstep(itype delta) noexcept {
        advance(-delta);
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    void engine<xtype, itype, output_mixin, output_previous, stream_mixin>::discard(itype delta) noexcept {
        advance(delta);
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    bool engine<xtype, itype, output_mixin, output_previous, stream_mixin>::wrapped() noexcept {
        if (stream_mixin::is_mcg) {
            // For MCGs, the low order two bits never change. In this
            // implementation, we keep them fixed at 3 to make this test
            // easier.
            return state_ == 3;
        }
        return state_ == 0;
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    template<typename... Args>
    void engine<xtype, itype, output_mixin, output_previous, stream_mixin>::seed(Args&&... args) noexcept {
        new (this) engine(std::forward<Args>(args)...);
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    itype engine<xtype, itype, output_mixin, output_previous, stream_mixin>::bump(itype state) noexcept {
        return state * multiplier<itype>() + increment();
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    itype engine<xtype, itype, output_mixin, output_previous, stream_mixin>::base_generate() noexcept {
        return state_ = bump(state_);
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    itype engine<xtype, itype, output_mixin, output_previous, stream_mixin>::base_generate0() noexcept {
        itype old_state = state_;
        state_ = bump(state_);
        return old_state;
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    itype engine<xtype, itype, output_mixin, output_previous, stream_mixin>::distance(itype newstate, itype mask) const noexcept {
        return distance(state_, newstate, multiplier<itype>(), increment(), mask);
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    itype engine<xtype, itype, output_mixin, output_previous, stream_mixin>::advance(
            itype state, itype delta, itype cur_mult, itype cur_plus) noexcept {
        // The method used here is based on Brown, "Random Number Generation
        // with Arbitrary Stride,", Transactions of the American Nuclear
        // Society (Nov. 1994).  The algorithm is very similar to fast
        // exponentiation.
        //
        // Even though delta is an unsigned integer, we can pass a
        // signed integer to go backwards, it just goes "the long way round".

        constexpr itype ZERO = 0U; // itype may be a non-trivial types, so
        constexpr itype ONE = 1U; // we define some ugly constants.
        itype acc_mult = 1;
        itype acc_plus = 0;
        while (delta > ZERO) {
            if (delta & ONE) {
                acc_mult *= cur_mult;
                acc_plus = acc_plus * cur_mult + cur_plus;
            }
            cur_plus = (cur_mult + ONE) * cur_plus;
            cur_mult *= cur_mult;
            delta >>= 1;
        }
        return acc_mult * state + acc_plus;
    }

    template<typename xtype, typename itype, typename output_mixin, bool output_previous, typename stream_mixin>
    itype engine<xtype, itype, output_mixin, output_previous, stream_mixin>::distance(
            itype cur_state, itype newstate, itype cur_mult, itype cur_plus, itype mask) noexcept {
        constexpr itype ONE = 1U; // itype could be weird, so use constant
        bool is_mcg = cur_plus == itype(0);
        itype the_bit = is_mcg ? itype(4U) : itype(1U);
        itype distance = 0U;
        while ((cur_state & mask) != (newstate & mask)) {
            if ((cur_state & the_bit) != (newstate & the_bit)) {
                cur_state = cur_state * cur_mult + cur_plus;
                distance |= the_bit;
            }
            assert((cur_state & the_bit) == (newstate & the_bit));
            the_bit <<= 1;
            cur_plus = (cur_mult + ONE) * cur_plus;
            cur_mult *= cur_mult;
        }
        return is_mcg ? distance >> 2 : distance;
    }
}
