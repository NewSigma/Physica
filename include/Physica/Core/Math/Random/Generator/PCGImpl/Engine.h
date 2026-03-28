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
#include <iosfwd>
#include <utility>
#include "OutputFn.h"

namespace Physica::Internal {
    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    class engine final : public Stream {
    public:
        using result_type = OutI;
        using state_type = WorkI;
    private:
        using OutFn = OutputFn<OutI, WorkI, Type>;
        using Stream::increment;

        WorkI state_;
    public:
        explicit engine(WorkI state = WorkI(0xCAFEF00DD15EA5E5ULL)) noexcept;
        template<typename sm = Stream>
        engine(WorkI state, typename sm::stream_state stream_seed) noexcept;
        explicit engine(auto& seedSeq) noexcept requires(!Stream::can_specify_stream);
        explicit engine(auto& seedSeq) noexcept requires(Stream::can_specify_stream);
        /* Operators */
        result_type operator()() noexcept;
        result_type operator()(result_type upper_bound) noexcept;

        template<typename xtype1, typename itype1, OutputFnType Type1, typename stream_mixin_lhs, typename stream_mixin_rhs>
        friend bool operator==(const engine<xtype1, itype1, Type1, stream_mixin_lhs>&,
                               const engine<xtype1, itype1, Type1, stream_mixin_rhs>&);

        template<typename xtype1, typename itype1, OutputFnType Type1, typename stream_mixin_lhs, typename stream_mixin_rhs>
        friend itype1 operator-(const engine<xtype1, itype1, Type1, stream_mixin_lhs>&,
                                const engine<xtype1, itype1, Type1, stream_mixin_rhs>&);

        template<typename CharT, typename Traits, typename xtype1, typename itype1, OutputFnType Type1, typename stream_mixin1>
        friend std::basic_ostream<CharT, Traits>&
        operator<<(std::basic_ostream<CharT, Traits>& out, const engine<xtype1, itype1, Type, stream_mixin1>&);

        template<typename CharT, typename Traits, typename xtype1, typename itype1, OutputFnType Type1, typename stream_mixin1>
        friend std::basic_istream<CharT, Traits>&
        operator>>(std::basic_istream<CharT, Traits>& in, engine<xtype1, itype1, Type, stream_mixin1>& rng);
        /* Operations */
        void advance(WorkI delta) noexcept;
        void backstep(WorkI delta) noexcept;
        void discard(WorkI delta) noexcept;

        template<typename... Args>
        void seed(Args&&... args) noexcept;
        /* Getters */
        [[nodiscard]] bool wrapped() const noexcept;
        /* Static members */
        constexpr static size_t period_pow2() noexcept { return sizeof(state_type) * 8 - 2 * Stream::is_mcg; }
        constexpr static result_type min() noexcept { return std::numeric_limits<result_type>::min(); }
        constexpr static result_type max() noexcept { return std::numeric_limits<result_type>::max(); }
    protected:
        WorkI bump(WorkI state) noexcept;
        WorkI base_generate0() noexcept;
        WorkI distance(WorkI newstate, WorkI mask = WorkI(~WorkI(0U))) const noexcept;
        /* Static members */
        static WorkI advance(WorkI state, WorkI delta, WorkI cur_mult, WorkI cur_plus) noexcept;
        static WorkI distance(WorkI cur_state, WorkI newstate, WorkI cur_mult, WorkI cur_plus, WorkI mask = ~WorkI(0U)) noexcept;
    };

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    engine<OutI, WorkI, Type, Stream>::engine(WorkI state) noexcept
            : state_(this->is_mcg ? (state | state_type(3U)) : bump(state + this->increment())) {}

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    template<typename sm>
    engine<OutI, WorkI, Type, Stream>::engine(WorkI state, typename sm::stream_state stream_seed) noexcept
            : Stream(stream_seed), state_(this->is_mcg ? state | state_type(3U) : bump(state + this->increment())) {}

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    engine<OutI, WorkI, Type, Stream>::engine(auto& seedSeq) noexcept requires(!Stream::can_specify_stream)
            : engine(seedSeq.template generate<WorkI>()) {}

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    engine<OutI, WorkI, Type, Stream>::engine(auto& seedSeq) noexcept requires(Stream::can_specify_stream) {
        std::array<WorkI, 2> seeddata{};
        seedSeq.generate(seeddata);
        seed(seeddata[1], seeddata[0]);
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    auto engine<OutI, WorkI, Type, Stream>::operator()() noexcept -> result_type {
        constexpr bool OutputPrev = sizeof(WorkI) <= 8;
        const WorkI old_state = state_;
        state_ = bump(state_);
        return OutFn::operator()(OutputPrev ? old_state : state_);
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    auto engine<OutI, WorkI, Type, Stream>::operator()(result_type upper_bound) noexcept -> result_type {
        return bounded_rand(*this, upper_bound);
    }

    template<typename CharT, typename Traits, class OutI, class WorkI, OutputFnType Type, class Stream>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& out, const engine<OutI, WorkI, Type, Stream>& rng) {
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

    template<typename CharT, typename Traits, class OutI, class WorkI, OutputFnType Type, class Stream>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& in, engine<OutI, WorkI, Type, Stream>& rng) {
        auto orig_flags = in.flags(std::ios_base::dec | std::ios_base::skipws);

        WorkI multiplier;
        WorkI increment;
        WorkI state;
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

    template<class OutI, class WorkI, OutputFnType Type, typename stream_mixin_lhs, typename stream_mixin_rhs>
    WorkI operator-(const engine<OutI, WorkI, Type, stream_mixin_lhs>& lhs,
                    const engine<OutI, WorkI, Type, stream_mixin_rhs>& rhs) {
        static_assert(std::same_as<stream_mixin_lhs, stream_mixin_rhs>, "Incomparable generators");
        if (lhs.increment() == rhs.increment())
            return rhs.distance(lhs.state_);

        constexpr WorkI ONE = 1U;
        WorkI lhs_diff = lhs.increment() + (lhs.multiplier() - ONE) * lhs.state_;
        WorkI rhs_diff = rhs.increment() + (rhs.multiplier() - ONE) * rhs.state_;
        if ((lhs_diff & WorkI(3U)) != (rhs_diff & WorkI(3U))) {
            rhs_diff = -rhs_diff;
        }
        return rhs.distance(rhs_diff, lhs_diff, rhs.multiplier(), WorkI(0U));
    }

    template<class OutI, class WorkI, OutputFnType Type, typename stream_mixin_lhs, typename stream_mixin_rhs>
    bool operator==(const engine<OutI, WorkI, Type, stream_mixin_lhs>& lhs,
                    const engine<OutI, WorkI, Type, stream_mixin_rhs>& rhs) {
        return (lhs.multiplier() == rhs.multiplier())
            && (lhs.increment() == rhs.increment())
            && (lhs.state_ == rhs.state_);
    }

    template<class OutI, class WorkI, OutputFnType Type, typename stream_mixin_lhs, typename stream_mixin_rhs>
    inline bool operator!=(const engine<OutI, WorkI, Type, stream_mixin_lhs>& lhs,
                           const engine<OutI, WorkI, Type, stream_mixin_rhs>& rhs) {
        return !operator==(lhs, rhs);
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    void engine<OutI, WorkI, Type, Stream>::advance(WorkI delta) noexcept {
        state_ = advance(state_, delta, this->multiplier(), this->increment());
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    void engine<OutI, WorkI, Type, Stream>::backstep(WorkI delta) noexcept {
        advance(-delta);
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    void engine<OutI, WorkI, Type, Stream>::discard(WorkI delta) noexcept {
        advance(delta);
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    template<typename... Args>
    void engine<OutI, WorkI, Type, Stream>::seed(Args&&... args) noexcept {
        new (this) engine(std::forward<Args>(args)...);
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    bool engine<OutI, WorkI, Type, Stream>::wrapped() const noexcept {
        // For MCGs, the low order two bits never change. In this
        // implementation, we keep them fixed at 3 to make this test
        // easier.
        if (Stream::is_mcg)
            return state_ == 3;
        return state_ == 0;
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    WorkI engine<OutI, WorkI, Type, Stream>::bump(WorkI state) noexcept {
        return state * OutFn::multiplier() + OutFn::increment();
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    WorkI engine<OutI, WorkI, Type, Stream>::base_generate0() noexcept {
        WorkI old_state = state_;
        state_ = bump(state_);
        return old_state;
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    WorkI engine<OutI, WorkI, Type, Stream>::distance(WorkI newstate, WorkI mask) const noexcept {
        return distance(state_, newstate, OutFn::multiplier(), OutFn::increment(), mask);
    }

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    WorkI engine<OutI, WorkI, Type, Stream>::advance(
            WorkI state, WorkI delta, WorkI cur_mult, WorkI cur_plus) noexcept {
        // The method used here is based on Brown, "Random Number Generation
        // with Arbitrary Stride,", Transactions of the American Nuclear
        // Society (Nov. 1994).  The algorithm is very similar to fast
        // exponentiation.
        //
        // Even though delta is an unsigned integer, we can pass a
        // signed integer to go backwards, it just goes "the long way round".

        constexpr WorkI ZERO = 0U; // WorkI may be a non-trivial types, so
        constexpr WorkI ONE = 1U; // we define some ugly constants.
        WorkI acc_mult = 1;
        WorkI acc_plus = 0;
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

    template<class OutI, class WorkI, OutputFnType Type, class Stream>
    WorkI engine<OutI, WorkI, Type, Stream>::distance(
            WorkI cur_state, WorkI newstate, WorkI cur_mult, WorkI cur_plus, WorkI mask) noexcept {
        constexpr WorkI ONE = 1U; // WorkI could be weird, so use constant
        bool is_mcg = cur_plus == WorkI(0);
        WorkI the_bit = is_mcg ? WorkI(4U) : WorkI(1U);
        WorkI distance = 0U;
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
