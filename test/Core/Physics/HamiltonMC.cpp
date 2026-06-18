/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/MC/HamiltonMC.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<>;

int main() {
    // Test that HMC does not get stuck on an empty model
    auto hmc = HamiltonMC<T>({1, 1});
    hmc.step<RandomSource>(EmptyForceModel<T, 1>{});
    return 0;
}
