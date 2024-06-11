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

namespace Physica::Core {
    template<unsigned int Dim>
    inline SpinlessFermion<Dim>::SpinlessFermion(uint64_t occupyBits_, int numSite_)
            : occupyBits(occupyBits_), numSite(numSite_) {
        assert((1 < numSite) && (static_cast<unsigned int>(numSite) < sizeof(occupyBits) * CHAR_BIT) && "[Error]: Invalid param");
    }

    template<unsigned int Dim>
    inline bool SpinlessFermion<Dim>::operator==(const This& other) const noexcept {
        return (occupyBits == other.occupyBits) && (numSite == other.numSite);
    }

    template<unsigned int Dim>
    inline SpinlessFermion<Dim> SpinlessFermion<Dim>::operator<<(int shift) const noexcept {
        assert(0 <= shift && shift < numSite && "[Error]: Invalid shift");
        const auto highBits = occupyBits << shift;
        const auto lowBits = occupyBits >> (numSite - shift);
        return SpinlessFermion((highBits | lowBits) & makeFullMask(), numSite);
    }

    template<unsigned int Dim>
    inline SpinlessFermion<Dim> SpinlessFermion<Dim>::operator>>(int shift) const noexcept {
        assert(0 <= shift && shift < numSite && "[Error]: Invalid shift");
        const auto highBits = occupyBits << (numSite - shift);
        const auto lowBits = occupyBits >> shift;
        return SpinlessFermion((highBits | lowBits) & makeFullMask(), numSite);
    }

    template<unsigned int Dim>
    inline SpinlessFermion<Dim> SpinlessFermion<Dim>::hop(unsigned char from, unsigned char to) const {
        const bool canHop = isOccupy(from) && !isOccupy(to);
        if (!canHop)
            return SpinlessFermion();
        const uint64_t fromMask = 1UL << from;
        const uint64_t toMask = 1UL << to;
        return SpinlessFermion((occupyBits ^ fromMask) | toMask, numSite);
    }

    template<unsigned int Dim>
    SpinlessFermion<Dim> SpinlessFermion<Dim>::transReduce(int period) const {
        assert(numSite % period == 0 && "[Error]: Invalid period");
        assert(0 < period && period <= numSite && "[Error]: Invalid period");
        if (period == numSite)
            return *this;
        const int numTrans = numSite / period;
        uint64_t result = occupyBits;
        This temp = *this;
        for (int i = 0; i < numTrans; ++i) {
            temp <<= period;
            result = std::min(result, temp.occupyBits);
        }
        return SpinlessFermion(result, numSite);
    }

    template<unsigned int Dim>
    int SpinlessFermion<Dim>::calcPeriod() const {
        This copy = *this;
        int i = 1;
        for (; i <= numSite; ++i) {
            copy <<= 1;
            if (copy == *this)
                break;
        }
        assert(i <= numSite && (numSite % i == 0) && "[Error]: Impossible, this is a bug");
        return i;
    }

    template<unsigned int Dim>
    inline void SpinlessFermion<Dim>::swap(SpinlessFermion& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(occupyBits, obj.occupyBits);
        std::swap(numSite, obj.numSite);
    }

    template<unsigned int Dim>
    inline bool SpinlessFermion<Dim>::isOccupy(unsigned char site) const noexcept {
        assert(site < sizeof(occupyBits) * CHAR_BIT && "[Error]: Invalid site");
        const uint64_t mask = 1UL << site;
        return (occupyBits & mask) != 0;
    }

    template<unsigned int Dim>
    inline bool SpinlessFermion<Dim>::isTransReducible(int period) const noexcept {
        return transReduce(period) != (*this);
    }

    template<unsigned int Dim>
    inline uint64_t SpinlessFermion<Dim>::makeFullMask() const noexcept {
        return (static_cast<uint64_t>(1) << numSite) - 1;
    }

    template<unsigned int Dim>
    inline uint64_t SpinlessFermion<Dim>::makeHighMask() const noexcept {
        return (static_cast<uint64_t>(1) << (numSite - 1));
    }
}
