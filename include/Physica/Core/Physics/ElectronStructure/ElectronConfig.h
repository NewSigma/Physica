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
#pragma once

#include "Physica/Utils/Container/Array/Array.h"

namespace Physica::Core::Physics {
    class ElectronConfig {
    public:
        enum OrbitState {
            NoOccupacy,
            SingleOccupacy,
            DoubleOccupacy
        };
    private:
        Utils::Array<OrbitState> states;
    public:
        ElectronConfig(size_t maxOrbitCount) : states(maxOrbitCount, NoOccupacy) { assert(maxOrbitCount > 0); }
        /* Setters */
        void setOrbitState(size_t orbitIndex, OrbitState state) { states[orbitIndex] = state; }
        /* Getters */
        [[nodiscard]] size_t getNumOrbit() const noexcept { return states.getLength(); }
        [[nodiscard]] OrbitState getOrbitState(size_t orbitIndex) const { return states[orbitIndex]; }
        [[nodiscard]] size_t getNumOccupiedOrbit() const noexcept;
        [[nodiscard]] size_t getOccupiedOrbitPos(size_t orbitIndex) const;
    };
}
