/*
 * Copyright 2024 WeiBo He.
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

#include <Physica/Core/Exception/NotImplementedException.h>

namespace Physica::Core {
    template<unsigned int Dim, unsigned int NumSite>
    inline SpinlessFermion<Dim, NumSite>::SpinlessFermion(IntType occupyBits_) : occupyBits(occupyBits_) {}

    template<unsigned int Dim, unsigned int NumSite>
    inline bool SpinlessFermion<Dim, NumSite>::operator==(const This& other) const noexcept {
        return occupyBits == other.occupyBits;
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline SpinlessFermion<Dim, NumSite> SpinlessFermion<Dim, NumSite>::operator<<(int shift) const noexcept {
        assert(0 <= shift && shift < int(NumSite) && "[Error]: Invalid shift");
        const auto highBits = occupyBits << shift;
        const auto lowBits = occupyBits >> (NumSite - shift);
        return SpinlessFermion((highBits | lowBits) & makeFullMask());
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline SpinlessFermion<Dim, NumSite> SpinlessFermion<Dim, NumSite>::operator>>(int shift) const noexcept {
        assert(0 <= shift && shift < int(NumSite) && "[Error]: Invalid shift");
        const auto highBits = occupyBits << (NumSite - shift);
        const auto lowBits = occupyBits >> shift;
        return SpinlessFermion((highBits | lowBits) & makeFullMask());
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline SpinlessFermion<Dim, NumSite> SpinlessFermion<Dim, NumSite>::hop(uint8_t from, uint8_t to) const {
        assert(from != to && "[Error]: Assuming different sites");
        const bool canHop = isOccupy(from) && !isOccupy(to);
        if (!canHop)
            return SpinlessFermion();
        const IntType fromMask = 1UL << from;
        const IntType toMask = 1UL << to;
        return SpinlessFermion((occupyBits ^ fromMask) | toMask);
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline int SpinlessFermion<Dim, NumSite>::hopSign(uint8_t from, uint8_t to) const {
        const int numElectron = countOnes(occupyBits >> (from + 1)) - countOnes(occupyBits >> (to + 1));
        const bool flag1 = from < to;
        const bool flag2 = numElectron % 2 == 0;
        return (flag1 == flag2) ? 1 : -1;
    }

    template<unsigned int Dim, unsigned int NumSite>
    SpinlessFermion<Dim, NumSite> SpinlessFermion<Dim, NumSite>::transReduce(int period) const {
        if constexpr (Dim != 1)
            throw NotImplementedException();
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
        return SpinlessFermion(result);
    }

    template<unsigned int Dim, unsigned int NumSite>
    int SpinlessFermion<Dim, NumSite>::lShiftSign() const {
        const This other = (*this) << 1;
        const bool noExchange = (*this) > other;
        if (noExchange)
            return 1;
        return (getNumParticle() % 2U != 0) ? 1 : -1;
    }

    template<unsigned int Dim, unsigned int NumSite>
    int SpinlessFermion<Dim, NumSite>::calcPeriod() const {
        if constexpr (Dim != 1)
            throw NotImplementedException();
        This copy = *this;
        unsigned int i = 1;
        for (; i <= NumSite; ++i) {
            copy <<= 1;
            if (copy == *this)
                break;
        }
        assert(i <= NumSite && (NumSite % i == 0) && "[Error]: Impossible, this is a bug");
        return i;
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline void SpinlessFermion<Dim, NumSite>::swap(SpinlessFermion& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(occupyBits, obj.occupyBits);
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline bool SpinlessFermion<Dim, NumSite>::isOccupy(uint8_t site) const noexcept {
        assert(site < sizeof(occupyBits) * CHAR_BIT && "[Error]: Invalid site");
        const IntType mask = 1UL << site;
        return (occupyBits & mask) != 0;
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline bool SpinlessFermion<Dim, NumSite>::isTransReducible(int period) const noexcept {
        return transReduce(period) != (*this);
    }

    template<unsigned int Dim, unsigned int NumSite>
    template<class RandomGenerator>
    SpinlessFermion<Dim, NumSite> SpinlessFermion<Dim, NumSite>::random_state(RandomGenerator& gen) {
        constexpr bool flag = NumSite < sizeof(typename RandomGenerator::result_type) * CHAR_BIT;
        static_assert(flag, "[Error]: The random generator cannot provide enough random bits");
        return gen() & makeFullMask();
    }

    template<unsigned int Dim, unsigned int NumSite>
    template<class RandomGenerator>
    SpinlessFermion<Dim, NumSite> SpinlessFermion<Dim, NumSite>::random_state(size_t numParticle, RandomGenerator& gen) {
        constexpr bool flag = NumSite < sizeof(typename RandomGenerator::result_type) * CHAR_BIT;
        static_assert(flag, "[Error]: The random generator cannot provide enough random bits");
        assert(numParticle <= NumSite);
        char bits[NumSite];
        memset(bits, 0, NumSite * sizeof(char));
        for (size_t i = 0; i < numParticle; ++i)
            bits[i] = 1;
        std::shuffle(bits, bits + NumSite, gen);

        IntType result = bits[0];
        for (unsigned int i = 1; i < NumSite; ++i) {
            result <<= 1U;
            result += bits[i];
        }
        return This(result);
    }
}
