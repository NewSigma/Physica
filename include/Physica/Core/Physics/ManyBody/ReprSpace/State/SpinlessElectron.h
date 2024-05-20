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

#include <iostream>
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

    class PHYSICA_API SpinlessElectron : public State<SpinlessElectron> {
        using This = SpinlessElectron;

        uint64_t occupyBits;
        int numSite;
    public:
        SpinlessElectron() = default;
        inline SpinlessElectron(uint64_t occupyBits_, int numSite_);
        SpinlessElectron(const SpinlessElectron&) = default;
        SpinlessElectron(SpinlessElectron&&) noexcept = default;
        ~SpinlessElectron() = default;
        /* Operators */
        This& operator=(This obj) noexcept { swap(obj); return *this; }
        [[nodiscard]] inline bool operator==(const This& other) const noexcept;
        [[nodiscard]] bool operator!=(const This& other) const noexcept { return !(*this == other); }
        [[nodiscard]] bool operator>(const This& other) const noexcept { return occupyBits > other.occupyBits; }
        [[nodiscard]] bool operator<(const This& other) const noexcept { return occupyBits < other.occupyBits; }
        [[nodiscard]] bool operator>=(const This& other) const noexcept { return occupyBits >= other.occupyBits; }
        [[nodiscard]] bool operator<=(const This& other) const noexcept { return occupyBits <= other.occupyBits; }
        [[nodiscard]] inline This operator<<(int shift) const noexcept;
        [[nodiscard]] inline This operator>>(int shift) const noexcept;
        void operator<<=(int shift) noexcept { (*this) = (*this) << shift; }
        void operator>>=(int shift) noexcept { (*this) = (*this) >> shift; }
        friend std::ostream& operator<<(std::ostream& os, SpinlessElectron e);
        /* Operations */
        [[nodiscard]] inline SpinlessElectron hop(unsigned char from, unsigned char to) const;
        [[nodiscard]] SpinlessElectron transReduce(int period = 1) const;
        [[nodiscard]] int calcPeriod() const;
        inline void swap(SpinlessElectron& __restrict obj) noexcept;
        /* Getters */
        [[nodiscard]] uint64_t getOccupyBits() const noexcept { return occupyBits; }
        [[nodiscard]] int getNumSite() const noexcept { return numSite; }
        [[nodiscard]] bool isVacuum() const noexcept { return occupyBits == 0; }
        [[nodiscard]] inline bool isOccupy(unsigned char site) const noexcept;
        [[nodiscard]] unsigned int getNumElectron() const noexcept { return countOnes(occupyBits); }
        [[nodiscard]] inline bool isTransReducible(int period = 1) const noexcept;
    private:
        [[nodiscard]] inline uint64_t makeFullMask() const noexcept;
        [[nodiscard]] inline uint64_t makeHighMask() const noexcept;
    };

    inline SpinlessElectron::SpinlessElectron(uint64_t occupyBits_, int numSite_)
            : occupyBits(occupyBits_), numSite(numSite_) {
        assert((1 < numSite) && (static_cast<unsigned int>(numSite) < sizeof(occupyBits) * CHAR_BIT) && "[Error]: Invalid param");
    }

    inline bool SpinlessElectron::operator==(const This& other) const noexcept {
        return (occupyBits == other.occupyBits) && (numSite == other.numSite);
    }

    inline SpinlessElectron SpinlessElectron::operator<<(int shift) const noexcept {
        assert(0 <= shift && shift < numSite && "[Error]: Invalid shift");
        const auto highBits = occupyBits << shift;
        const auto lowBits = occupyBits >> (numSite - shift);
        return SpinlessElectron((highBits | lowBits) & makeFullMask(), numSite);
    }

    inline SpinlessElectron SpinlessElectron::operator>>(int shift) const noexcept {
        assert(0 <= shift && shift < numSite && "[Error]: Invalid shift");
        const auto highBits = occupyBits << (numSite - shift);
        const auto lowBits = occupyBits >> shift;
        return SpinlessElectron((highBits | lowBits) & makeFullMask(), numSite);
    }

    inline SpinlessElectron SpinlessElectron::hop(unsigned char from, unsigned char to) const {
        const bool canHop = isOccupy(from) && !isOccupy(to);
        if (!canHop)
            return SpinlessElectron();
        const uint64_t fromMask = 1UL << from;
        const uint64_t toMask = 1UL << to;
        return SpinlessElectron((occupyBits ^ fromMask) | toMask, numSite);
    }

    inline void SpinlessElectron::swap(SpinlessElectron& __restrict obj) noexcept {
        assert(this != &obj && "[Error]: Self swap is likely a bug");
        std::swap(occupyBits, obj.occupyBits);
        std::swap(numSite, obj.numSite);
    }

    inline bool SpinlessElectron::isOccupy(unsigned char site) const noexcept {
        assert(site < sizeof(occupyBits) * CHAR_BIT && "[Error]: Invalid site");
        const uint64_t mask = 1UL << site;
        return (occupyBits & mask) != 0;
    }

    inline bool SpinlessElectron::isTransReducible(int period) const noexcept {
        return transReduce(period) != (*this);
    }

    inline uint64_t SpinlessElectron::makeFullMask() const noexcept {
        return (static_cast<uint64_t>(1) << numSite) - 1;
    }

    inline uint64_t SpinlessElectron::makeHighMask() const noexcept {
        return (static_cast<uint64_t>(1) << (numSite - 1));
    }

    inline std::ostream& operator<<(std::ostream& os, SpinlessElectron e) {
        auto mask = e.makeHighMask();
        for (int i = 0; i < e.getNumSite(); ++i) {
            const bool flag = (e.getOccupyBits() & mask) == 0;
            os << (flag ? '0' : '1');
            mask >>= 1;
        }
        return os;
    }
}
