/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Physics/ManyBody/ReprSpace/State/SpinState.h"
#include "Test.h"

using namespace Physica;

int main() {
    const SpinState<1, 6> psi(0b011010);
    const auto psi1 = psi << 2;
    expect(psi1.getOccupyBits() == 0b101001);
    expect((psi1 >> 2) == psi);
    expect(!psi.isTransReducible(2));
    return 0;
}
