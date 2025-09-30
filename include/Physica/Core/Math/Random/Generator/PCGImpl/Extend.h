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

#include <array>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include "OutputFn.h"

namespace Physica::Internal {
    template<typename baseclass>
    struct inside_out : private baseclass {
        inside_out() = delete;

        using result_type = typename baseclass::result_type;
        using state_type = typename baseclass::state_type;
        static_assert(sizeof(result_type) == sizeof(state_type),
                      "Require a RNG whose output function is a permutation");

        static bool external_step(result_type& randval, size_t i) {
            state_type state = baseclass::unoutput(randval);
            state = state * baseclass::multiplier() + baseclass::increment()
                  + state_type(i * 2);
            result_type result = baseclass::output(state);
            randval = result;
            state_type zero =
                    baseclass::is_mcg ? state & state_type(3U) : state_type(0U);
            return result == zero;
        }

        static bool external_advance(result_type& randval, size_t i, result_type delta, bool forwards = true) {
            state_type state = baseclass::unoutput(randval);
            state_type mult = baseclass::multiplier();
            state_type inc = baseclass::increment() + state_type(i * 2);
            state_type zero =
                    baseclass::is_mcg ? state & state_type(3U) : state_type(0U);
            state_type dist_to_zero = baseclass::distance(state, zero, mult, inc);
            bool crosses_zero =
                    forwards ? dist_to_zero <= delta
                             : (-dist_to_zero) <= delta;
            if (!forwards)
                delta = -delta;
            state = baseclass::advance(state, delta, mult, inc);
            randval = baseclass::output(state);
            return crosses_zero;
        }
    };

    template<uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD = true>
    class extended : public baseclass {
        using This = extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>;
    public:
        using state_type = typename baseclass::state_type;
        using result_type = typename baseclass::result_type;
        using insideout = inside_out<extvalclass>;
    private:
        constexpr static uint8_t rtypebits = sizeof(result_type) * 8;
        constexpr static uint8_t stypebits = sizeof(state_type) * 8;

        constexpr static uint8_t tick_limit_pow2 = 64U;

        constexpr static size_t table_size = 1UL << table_pow2;
        constexpr static size_t table_shift = stypebits - table_pow2;
        constexpr static state_type table_mask =
                (state_type(1U) << table_pow2) - state_type(1U);

        constexpr static bool may_tick =
                (advance_pow2 < stypebits) && (advance_pow2 < tick_limit_pow2);
        constexpr static size_t tick_shift = stypebits - advance_pow2;
        constexpr static state_type tick_mask =
                may_tick ? state_type(
                                   (uint64_t(1) << (advance_pow2 * may_tick)) - 1)
                         // ^-- stupidity to appease GCC warnings
                         : ~state_type(0U);

        constexpr static bool may_tock = stypebits < tick_limit_pow2;

        std::array<result_type, table_size> data_;

        void advance_table();

        void advance_table(state_type delta, bool isForwards = true);

        result_type& get_extended_value() {
            state_type state = this->state_;
            if (KDD && baseclass::is_mcg) {
                // The low order bits of an MCG are constant, so drop them.
                state >>= 2;
            }
            size_t index = KDD ? state & table_mask
                               : state >> table_shift;

            if (may_tick) {
                bool tick = KDD ? (state & tick_mask) == state_type(0U)
                                : (state >> tick_shift) == state_type(0U);
                if (tick)
                    advance_table();
            }
            if (may_tock) {
                bool tock = state == state_type(0U);
                if (tock)
                    advance_table();
            }
            return data_[index];
        }
    public:
        template<typename SeedSeq>
        extended(SeedSeq& seedSeq) requires(!std::is_convertible<SeedSeq, result_type>::value && !std::is_convertible<SeedSeq, extended>::value)
                : baseclass(seedSeq) {
            seedSeq.generate(data_);
        }

        extended(const result_type* data)
                : baseclass() {
            datainit(data);
        }

        extended(const result_type* data, state_type seed)
                : baseclass(seed) {
            datainit(data);
        }

        // This function may or may not exist.  It thus has to be a template
        // to use SFINAE; users don't have to worry about its template-ness.

        template<typename bc = baseclass>
        extended(const result_type* data, state_type seed, typename bc::stream_state stream_seed)
                : baseclass(seed, stream_seed) {
            datainit(data);
        }

        extended()
                : baseclass() {
            selfinit();
        }

        extended(state_type seed)
                : baseclass(seed) {
            selfinit();
        }

        // This function may or may not exist.  It thus has to be a template
        // to use SFINAE; users don't have to worry about its template-ness.

        template<typename bc = baseclass>
        extended(state_type seed, typename bc::stream_state stream_seed)
                : baseclass(seed, stream_seed) {
            selfinit();
        }
        /* Operators */
        result_type operator()() {
            result_type rhs = get_extended_value();
            result_type lhs = this->baseclass::operator()();
            return lhs ^ rhs;
        }

        result_type operator()(result_type upper_bound) {
            return bounded_rand(*this, upper_bound);
        }

        bool operator==(const This& other) const;

        template<typename CharT, typename Traits, uint8_t table_pow2_, uint8_t advance_pow2_, typename baseclass_, typename extvalclass_, bool KDD_>
        friend std::basic_ostream<CharT, Traits>&
        operator<<(std::basic_ostream<CharT, Traits>& out, const extended<table_pow2_, advance_pow2_, baseclass_, extvalclass_, KDD_>&);

        template<typename CharT, typename Traits, uint8_t table_pow2_, uint8_t advance_pow2_, typename baseclass_, typename extvalclass_, bool KDD_>
        friend std::basic_istream<CharT, Traits>&
        operator>>(std::basic_istream<CharT, Traits>& in, extended<table_pow2_, advance_pow2_, baseclass_, extvalclass_, KDD_>&);
        /* Operations */
        template<typename... Args>
        void seed(Args&&... args) {
            new (this) extended(std::forward<Args>(args)...);
        }

        void set(result_type wanted) {
            result_type& rhs = get_extended_value();
            result_type lhs = this->baseclass::operator()();
            rhs = lhs ^ wanted;
        }

        void advance(state_type distance, bool forwards = true);

        void backstep(state_type distance) {
            advance(distance, false);
        }
        /* Getters */
        constexpr static size_t period_pow2() { return baseclass::period_pow2() + table_size * extvalclass::period_pow2(); }
    private:
        void selfinit();
        void datainit(const result_type* data);
    };

    template<uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD>
    void extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>::datainit(const result_type* data) {
        for (size_t i = 0; i < table_size; ++i)
            data_[i] = data[i];
    }

    template<uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD>
    void extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>::selfinit() {
        // We need to fill the extended table with something, and we have
        // very little provided data, so we use the base generator to
        // produce values.  Although not ideal (use a seed sequence, folks!),
        // unexpected correlations are mitigated by
        //      - using XOR differences rather than the number directly
        //      - the way the table is accessed, its values *won't* be accessed
        //        in the same order the were written.
        //      - any strange correlations would only be apparent if we
        //        were to backstep the generator so that the base generator
        //        was generating the same values again
        result_type lhs = baseclass::operator()();
        result_type rhs = baseclass::operator()();
        result_type xdiff = lhs - rhs;
        for (size_t i = 0; i < table_size; ++i) {
            data_[i] = baseclass::operator()() ^ xdiff;
        }
    }

    template<uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD>
    bool extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>::operator==(const This& other) const {
        auto& base_lhs = static_cast<const baseclass&>(*this);
        auto& base_rhs = static_cast<const baseclass&>(other);
        return base_lhs == base_rhs
            && std::equal(data_.begin(), data_.end(), other.data_.begin());
    }

    template<uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD>
    inline bool operator!=(const extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>& lhs,
                           const extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>& rhs) {
        return !operator==(lhs, rhs);
    }

    template<typename CharT, typename Traits, uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& out,
               const extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>& rng) {
        auto orig_flags = out.flags(std::ios_base::dec | std::ios_base::left);
        auto space = out.widen(' ');
        auto orig_fill = out.fill();

        out << rng.multiplier() << space
            << rng.increment() << space
            << rng.state_;

        for (const auto& datum : rng.data_)
            out << space << datum;

        out.flags(orig_flags);
        out.fill(orig_fill);
        return out;
    }

    template<typename CharT, typename Traits, uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& in, extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>& rng) {
        extended<table_pow2, advance_pow2, baseclass, extvalclass> new_rng;
        auto& base_rng = static_cast<baseclass&>(new_rng);
        in >> base_rng;

        if (in.fail())
            return in;

        auto orig_flags = in.flags(std::ios_base::dec | std::ios_base::skipws);

        for (auto& datum : new_rng.data_) {
            in >> datum;
            if (in.fail())
                break;
        }
        rng = new_rng;
        in.flags(orig_flags);
        return in;
    }

    template<uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD>
    void extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>::advance_table() {
        bool carry = false;
        for (size_t i = 0; i < table_size; ++i) {
            if (carry) {
                carry = insideout::external_step(data_[i], i + 1);
            }
            bool carry2 = insideout::external_step(data_[i], i + 1);
            carry = carry || carry2;
        }
    }

    template<uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD>
    void extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>::advance_table(
            state_type delta, bool isForwards) {
        using base_state_t = typename baseclass::state_type;
        using ext_state_t = typename extvalclass::state_type;
        constexpr uint8_t basebits = sizeof(base_state_t) * 8;
        constexpr uint8_t extbits = sizeof(ext_state_t) * 8;
        static_assert(basebits <= extbits || advance_pow2 > 0,
                      "Current implementation might overflow its carry");

        base_state_t carry = 0;
        for (size_t i = 0; i < table_size; ++i) {
            base_state_t total_delta = carry + delta;
            ext_state_t trunc_delta = ext_state_t(total_delta);
            if (basebits > extbits) {
                carry = total_delta >> extbits;
            }
            else {
                carry = 0;
            }
            carry += insideout::external_advance(data_[i], i + 1, trunc_delta, isForwards);
        }
    }

    template<uint8_t table_pow2, uint8_t advance_pow2, typename baseclass, typename extvalclass, bool KDD>
    void extended<table_pow2, advance_pow2, baseclass, extvalclass, KDD>::advance(
            state_type distance, bool forwards) {
        static_assert(KDD,
                      "Efficient advance is too hard for non-KDD extension. "
                      "For a weak advance, cast to base class");
        state_type zero =
                baseclass::is_mcg ? this->state_ & state_type(3U) : state_type(0U);
        if (may_tick) {
            state_type ticks = distance >> (advance_pow2 * may_tick);
            // ^-- stupidity to appease GCC
            // warnings
            state_type adv_mask =
                    baseclass::is_mcg ? tick_mask << 2 : tick_mask;
            state_type next_advance_distance = this->distance(zero, adv_mask);
            if (!forwards)
                next_advance_distance = (-next_advance_distance) & tick_mask;
            if (next_advance_distance < (distance & tick_mask)) {
                ++ticks;
            }
            if (ticks)
                advance_table(ticks, forwards);
        }
        if (forwards) {
            if (may_tock && this->distance(zero) <= distance)
                advance_table();
            baseclass::advance(distance);
        }
        else {
            if (may_tock && -(this->distance(zero)) <= distance)
                advance_table(state_type(1U), false);
            baseclass::advance(-distance);
        }
    }
}
