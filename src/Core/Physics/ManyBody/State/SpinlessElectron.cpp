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
#include "Physica/Core/Physics/ManyBody/ReprSpace/State/SpinlessElectron.h"

namespace Physica::Core {
    SpinlessElectron::SpinlessElectron(uint64_t occupyBits_, int numSite_)
            : occupyBits(occupyBits_), numSite(numSite_) {
        assert((1 < numSite) && (static_cast<unsigned int>(numSite) < sizeof(occupyBits) * CHAR_BIT) && "[Error]: Invalid param");
    }

    SpinlessElectron SpinlessElectron::operator<<(int shift) const {
        assert(0 <= shift && shift < numSite && "[Error]: Invalid shift");
        const auto highMask = makeHighMask();
        auto bits = occupyBits;
        for (int i = 0; i < shift; ++i) {
            const bool highBit = (bits & highMask) != 0;
            bits = (bits << 1) | highBit;
        }
        return SpinlessElectron(bits & makeFullMask(), numSite);
    }

    SpinlessElectron SpinlessElectron::operator>>(int shift) const {
        assert(0 <= shift && shift < numSite && "[Error]: Invalid shift");
        const auto highMask = makeHighMask();
        auto bits = occupyBits;
        for (int i = 0; i < shift; ++i) {
            const bool lowBit = (bits & static_cast<uint64_t>(1)) != 0;
            bits >>= 1;
            bits |= lowBit ? highMask : 0;
        }
        return SpinlessElectron(bits, numSite);
    }

    SpinlessElectron SpinlessElectron::transReduce(int period) const {
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
        return SpinlessElectron(result, numSite);
    }

    int SpinlessElectron::calcPeriod() const {
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
}
