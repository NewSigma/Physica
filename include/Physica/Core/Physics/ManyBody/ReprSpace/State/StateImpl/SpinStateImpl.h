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

#include "Physica/Core/Exception/NoImplException.h"
#include "../SpinState.h"

namespace Physica {
    template<int Dim, int NumSite>
    inline SpinState<Dim, NumSite>::SpinState(IntType occupyBits_) : occupyBits(occupyBits_) {}

    template<int Dim, int NumSite>
    inline bool SpinState<Dim, NumSite>::operator==(const This& other) const noexcept {
        return occupyBits == other.occupyBits;
    }

    template<int Dim, int NumSite>
    inline SpinState<Dim, NumSite> SpinState<Dim, NumSite>::operator<<(int shift) const noexcept {
        assert(0 <= shift && shift < int(NumSite) && "[Error]: Invalid shift");
        const auto highBits = occupyBits << shift;
        const auto lowBits = occupyBits >> (NumSite - shift);
        return SpinState((highBits | lowBits) & makeFullMask());
    }

    template<int Dim, int NumSite>
    inline SpinState<Dim, NumSite> SpinState<Dim, NumSite>::operator>>(int shift) const noexcept {
        assert(0 <= shift && shift < int(NumSite) && "[Error]: Invalid shift");
        const auto highBits = occupyBits << (NumSite - shift);
        const auto lowBits = occupyBits >> shift;
        return SpinState((highBits | lowBits) & makeFullMask());
    }

    template<int Dim, int NumSite>
    inline auto SpinState<Dim, NumSite>::operator&(const This& psi) const noexcept -> This {
        return This(occupyBits & psi.occupyBits);
    }

    template<int Dim, int NumSite>
    inline auto SpinState<Dim, NumSite>::operator|(const This& psi) const noexcept -> This {
        return This(occupyBits | psi.occupyBits);
    }

    template<int Dim, int NumSite>
    inline auto SpinState<Dim, NumSite>::operator^(const This& psi) const noexcept -> This {
        return This(occupyBits ^ psi.occupyBits);
    }

    template<int Dim, int NumSite>
    inline auto SpinState<Dim, NumSite>::operator~() const noexcept -> This {
        return This(~occupyBits);
    }

    template<int Dim, int NumSite>
    auto SpinState<Dim, NumSite>::flip(int8_t site) const -> This {
        assert(site < NumSite);
        return This(occupyBits ^ (IntType(1) << site));
    }

    template<int Dim, int NumSite>
    inline auto SpinState<Dim, NumSite>::hop(int8_t from, int8_t to) const -> This {
        assert(from != to && "[Error]: Assuming different sites");
        const bool canHop = isOccupy(from) && !isOccupy(to);
        if (!canHop)
            return SpinState();
        const IntType fromMask = 1UL << from;
        const IntType toMask = 1UL << to;
        return SpinState((occupyBits ^ fromMask) | toMask);
    }

    template<int Dim, int NumSite>
    inline int SpinState<Dim, NumSite>::hopSign(int8_t from, int8_t to) const {
        const int numElectron = std::popcount(occupyBits >> (from + 1)) - std::popcount(occupyBits >> (to + 1));
        const bool flag1 = from < to;
        const bool flag2 = numElectron % 2 == 0;
        return (flag1 == flag2) ? 1 : -1;
    }

    template<int Dim, int NumSite>
    auto SpinState<Dim, NumSite>::transReduce(int period) const -> This {
        if constexpr (Dim != 1)
            noImpl(__func__);
        assert(NumSite % period == 0 && "[Error]: Invalid period");
        assert(0 < period && period <= int(NumSite) && "[Error]: Invalid period");
        if (period == NumSite)
            return *this;
        const int numTrans = NumSite / period;
        IntType result = occupyBits;
        This temp = *this;
        for (int i = 0; i < numTrans; ++i) {
            temp <<= period;
            result = std::min(result, temp.occupyBits);
        }
        return SpinState(result);
    }

    template<int Dim, int NumSite>
    int SpinState<Dim, NumSite>::lShiftSign() const {
        const This other = (*this) << 1;
        const bool noExchange = (*this) > other;
        if (noExchange)
            return 1;
        return (getNumParticle() % 2U != 0) ? 1 : -1;
    }

    template<int Dim, int NumSite>
    int SpinState<Dim, NumSite>::calcPeriod() const {
        if constexpr (Dim != 1)
            noImpl(__func__);
        This copy = *this;
        int i = 1;
        for (; i <= NumSite; ++i) {
            copy <<= 1;
            if (copy == *this)
                break;
        }
        assert(i <= NumSite && (NumSite % i == 0) && "[Error]: Impossible, this is a bug");
        return i;
    }

    template<int Dim, int NumSite>
    inline void SpinState<Dim, NumSite>::swap(SpinState& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(occupyBits, obj.occupyBits);
    }

    template<int Dim, int NumSite>
    inline bool SpinState<Dim, NumSite>::isOccupy(int8_t site) const noexcept {
        assert((0 <= site) && (static_cast<uint8_t>(site) < sizeof(occupyBits) * CHAR_BIT) && "[Error]: Invalid site");
        const IntType mask = 1UL << site;
        return (occupyBits & mask) != 0;
    }

    template<int Dim, int NumSite>
    inline bool SpinState<Dim, NumSite>::isTransReducible(int period) const noexcept {
        return transReduce(period) != (*this);
    }

    template<int Dim, int NumSite>
    inline bool SpinState<Dim, NumSite>::isSpinUp(int8_t site) const noexcept {
        return isOccupy(site);
    }

    template<int Dim, int NumSite>
    inline bool SpinState<Dim, NumSite>::isSpinDown(int8_t site) const noexcept {
        return !isSpinUp(site);
    }

    template<int Dim, int NumSite>
    template<RNG R>
    SpinState<Dim, NumSite> SpinState<Dim, NumSite>::random_state() {
        constexpr bool flag = NumSite < sizeof(typename R::result_type) * CHAR_BIT;
        static_assert(flag, "[Error]: The random generator cannot provide enough random bits");
        return R::getInstance()() & makeFullMask();
    }

    template<int Dim, int NumSite>
    template<RNG R>
    SpinState<Dim, NumSite> SpinState<Dim, NumSite>::random_state(size_t numParticle) {
        constexpr bool flag = NumSite < sizeof(typename R::result_type) * CHAR_BIT;
        static_assert(flag, "[Error]: The random generator cannot provide enough random bits");
        assert(numParticle <= NumSite);
        char bits[NumSite];
        memset(bits, 0, NumSite * sizeof(char));
        for (size_t i = 0; i < numParticle; ++i)
            bits[i] = 1;
        std::shuffle(bits, bits + NumSite, R::getInstance());

        IntType result = bits[0];
        for (int i = 1; i < NumSite; ++i) {
            result <<= 1U;
            result += bits[i];
        }
        return This(result);
    }
}
