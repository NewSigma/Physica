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
    inline SpinlessFermion<Dim, NumSite>::SpinlessFermion(uint64_t occupyBits_) : occupyBits(occupyBits_) {}

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
    inline SpinlessFermion<Dim, NumSite> SpinlessFermion<Dim, NumSite>::hop(unsigned char from, unsigned char to) const {
        assert(from != to && "[Error]: Assuming different sites");
        const bool canHop = isOccupy(from) && !isOccupy(to);
        if (!canHop)
            return SpinlessFermion();
        const uint64_t fromMask = 1UL << from;
        const uint64_t toMask = 1UL << to;
        return SpinlessFermion((occupyBits ^ fromMask) | toMask);
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline SpinlessFermion<Dim, NumSite> SpinlessFermion<Dim, NumSite>::hop(IndexType dims, IndexType from, IndexType to) const {
        return hop(IndexType::toIndex1D(dims, from), IndexType::toIndex1D(dims, to));
    }

    template<unsigned int Dim, unsigned int NumSite>
    int SpinlessFermion<Dim, NumSite>::hopSign(unsigned char from, unsigned char to) const {
        const bool flag = from < to;
        const int sign1 = flag ? 1 : -1;
        if (!flag)
            std::swap(from, to);
        const unsigned int numElectron = countOnes(occupyBits >> (from + 1)) - countOnes(occupyBits >> (to + 1));
        const int sign2 = (numElectron % 2U == 0U) ? 1 : -1;
        return sign1 * sign2;
    }

    template<unsigned int Dim, unsigned int NumSite>
    int SpinlessFermion<Dim, NumSite>::hopSign(IndexType dims, IndexType from, IndexType to) const {
        return hopSign(IndexType::toIndex1D(dims, from), IndexType::toIndex1D(dims, to));
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
        uint64_t result = occupyBits;
        This temp = *this;
        for (int i = 0; i < numTrans; ++i) {
            temp <<= period;
            result = std::min(result, temp.occupyBits);
        }
        return SpinlessFermion(result);
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
    inline bool SpinlessFermion<Dim, NumSite>::isOccupy(unsigned char site) const noexcept {
        assert(site < sizeof(occupyBits) * CHAR_BIT && "[Error]: Invalid site");
        const uint64_t mask = 1UL << site;
        return (occupyBits & mask) != 0;
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline bool SpinlessFermion<Dim, NumSite>::isTransReducible(int period) const noexcept {
        return transReduce(period) != (*this);
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline uint64_t SpinlessFermion<Dim, NumSite>::makeFullMask() const noexcept {
        return (static_cast<uint64_t>(1) << NumSite) - 1;
    }

    template<unsigned int Dim, unsigned int NumSite>
    inline uint64_t SpinlessFermion<Dim, NumSite>::makeHighMask() const noexcept {
        return (static_cast<uint64_t>(1) << (NumSite - 1));
    }
}
