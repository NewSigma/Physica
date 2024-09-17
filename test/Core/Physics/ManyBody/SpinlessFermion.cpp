/*
 * Copyright 2024 Weibo He.
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
#include <iostream>
#include <Physica/Core/Physics/ManyBody/ReprSpace/State/SpinlessFermion.h>

using namespace Physica::Core;

int main() {
    const SpinlessFermion<1, 6> psi(0b011010);
    const auto psi1 = psi << 2;
    if (psi1.getOccupyBits() != 0b101001)
        return 1;
    if ((psi1 >> 2) != psi)
        return 1;
    if (psi.isTransReducible(2))
        return 1;
    return 0;
}
