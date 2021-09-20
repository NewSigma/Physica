/*
 * Copyright 2021 WeiBo He.
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
#include "Physica/Core/Physics/ElectronStructure/ElectronConfig.h"

namespace Physica::Core::Physics {
    size_t ElectronConfig::getNumOccupiedOrbit() const noexcept {
        size_t result = 0;
        for (size_t i = 0; i < states.getLength(); ++i)
            result += states[i] != NoOccupacy;
        return result;
    }

    size_t ElectronConfig::getOccupiedOrbitPos(size_t orbitIndex) const {
        assert(orbitIndex < getNumOccupiedOrbit());
        size_t i = 0;
        size_t scannedOrbit = 0;
        for (; i < states.getLength() && scannedOrbit <= orbitIndex; ++i)
            scannedOrbit += states[i] != NoOccupacy;
        assert(scannedOrbit == orbitIndex + 1);
        return i - 1;
    }
}
