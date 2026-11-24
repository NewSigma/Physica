/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2024 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "Physica/Macro.h"
#include "Physica/Core/Utils/Container/Array.h"

namespace Physica {
    class PHYSICA_API ElectronConfig {
    public:
        enum OrbitState {
            NoOccupacy,
            SingleOccupacy,
            DoubleOccupacy
        };
    private:
        Array<OrbitState> states;
    public:
        ElectronConfig(size_t maxOrbitCount) : states(maxOrbitCount, NoOccupacy) { assert(maxOrbitCount > 0); }
        /* Setters */
        void setOrbitState(size_t orbitIndex, OrbitState state) { states[orbitIndex] = state; }
        /* Getters */
        [[nodiscard]] size_t getNumOrbit() const noexcept { return states.getLength(); }
        [[nodiscard]] OrbitState getOrbitState(size_t orbitIndex) const { return states[orbitIndex]; }
        [[nodiscard]] size_t getNumOccupiedOrbit() const noexcept;
        [[nodiscard]] size_t getOccupiedOrbitPos(size_t orbitIndex) const;
        [[nodiscard]] size_t getElectronCount() const;
        [[nodiscard]] size_t getUpSpinElectronCount() const noexcept { return getNumOccupiedOrbit(); }
        [[nodiscard]] size_t getDownSpinElectronCount() const;
    };
}
