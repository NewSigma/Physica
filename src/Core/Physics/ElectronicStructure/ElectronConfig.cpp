/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Core/Physics/ElectronicStructure/ElectronConfig.h"

using namespace Physica;

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

size_t ElectronConfig::getElectronCount() const {
    size_t result = 0;
    for (size_t i = 0; i < states.getLength(); ++i) {
        switch (states[i]) {
        case NoOccupacy:
            break;
        case SingleOccupacy:
            result += 1;
            break;
        case DoubleOccupacy:
            result += 2;
            break;
        }
    }
    return result;
}

size_t ElectronConfig::getDownSpinElectronCount() const {
    size_t result = 0;
    for (size_t i = 0; i < states.getLength(); ++i)
        result += states[i] == DoubleOccupacy;
    return result;
}
