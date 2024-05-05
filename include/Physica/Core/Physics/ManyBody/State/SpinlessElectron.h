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

#include <iosfwd>
#include "State.h"

namespace Physica::Core {
    class SpinlessElectron;

    namespace Internal {
        template<>
        class Traits<SpinlessElectron> {
        public:
            constexpr static unsigned int SiteDOF = 2;
        };
    }

    class SpinlessElectron : public State<SpinlessElectron> {
        uint64_t occupyBits;
    public:
        SpinlessElectron() : occupyBits(0) {}
        SpinlessElectron(uint64_t occupyBits_) : occupyBits(occupyBits_) {}
        SpinlessElectron(const SpinlessElectron&) = default;
        SpinlessElectron(SpinlessElectron&&) noexcept = default;
        ~SpinlessElectron() = default;
        /* Operators */
        SpinlessElectron& operator=(SpinlessElectron obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] bool operator==(SpinlessElectron other) const noexcept { return occupyBits == other.occupyBits; }
        /* Operations */
        [[nodiscard]] inline SpinlessElectron hop(unsigned char from, unsigned char to) const;
        [[nodiscard]] inline SpinlessElectron transReduce(unsigned char numSite) const;
        inline void swap(SpinlessElectron& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] uint64_t getOccupyBits() const noexcept { return occupyBits; }
        [[nodiscard]] bool isVacuum() const noexcept { return occupyBits == 0; }
        [[nodiscard]] inline bool isOccupy(unsigned char site) const noexcept;
        [[nodiscard]] unsigned int getNumElectron() const noexcept { return countOnes(occupyBits); }
        [[nodiscard]] inline bool isTransIrreducible(unsigned char numSite) const noexcept;
    };

    inline SpinlessElectron SpinlessElectron::hop(unsigned char from, unsigned char to) const {
        const bool canHop = isOccupy(from) && !isOccupy(to);
        if (!canHop)
            return SpinlessElectron();
        const uint64_t fromMask = 1UL << from;
        const uint64_t toMask = 1UL << to;
        return SpinlessElectron((occupyBits ^ fromMask) | toMask);
    }

    inline SpinlessElectron SpinlessElectron::transReduce(unsigned char numSite) const {
        assert((1 < numSite) && (numSite < sizeof(occupyBits) * CHAR_BIT) && "[Error]: Invalid param");
        uint64_t result = occupyBits;
        const uint64_t fullMask = (static_cast<uint64_t>(1) << numSite) - 1;
        if (isVacuum() || result == fullMask)
            return result;

        const uint64_t highMask = (static_cast<uint64_t>(1) << (numSite - 1));
        while ((result & highMask) != 0) {
            result = (result << 1U) + 1U;
            result &= fullMask;
        }

        while (result % 2U == 0U)
            result >>= 1U;
        return result;
    }

    inline void SpinlessElectron::swap(SpinlessElectron& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(occupyBits, obj.occupyBits);
    }

    inline bool SpinlessElectron::isOccupy(unsigned char site) const noexcept {
        assert(site < sizeof(occupyBits) * CHAR_BIT && "[Error]: Invalid site");
        const uint64_t mask = 1UL << site;
        return (occupyBits & mask) != 0;
    }

    inline bool SpinlessElectron::isTransIrreducible(unsigned char numSite) const noexcept {
        assert((1 < numSite) && (numSite < sizeof(occupyBits) * CHAR_BIT) && "[Error]: Invalid param");
        const bool isOdd = occupyBits % 2U == 1U;
        const bool isLowerHalf = static_cast<uint64_t>(1) << (numSite - 1);
        const bool isFullOccupy = (static_cast<uint64_t>(1) << numSite) - 1 == occupyBits;
        return (isOdd && isLowerHalf) || isFullOccupy;
    }

    inline std::ostream& operator<<(std::ostream& os, SpinlessElectron e) {
        return os << e.getOccupyBits();
    }
}
