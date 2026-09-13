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
#include "Test.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<>;

namespace {
    void empty() {
        // Test that HMC does not get stuck on an empty model
        auto hmc = HamiltonMC<T>({1, 1});
        hmc.step<RandomSource>(EmptyForceModel<T, 1>{});
    }

    void initialize() {
        // Test that temperature is not fixed
        using KineticModel = OpenModel<T, 1, 1>;
        constexpr size_t numSample = 1024;
        auto hmc = HamiltonMC<T>({1, 1});
        auto& root = hmc.getRoot();
        T mean = 0;
        T var = 0;
        for (size_t i = 0; i < numSample; ++i) {
            root.template initMomentum<KineticModel, RandomSource>();
            var.toNextVariance(mean, i, root.template calcTemperature<KineticModel>());
        }
        expect<RandomSource>(scalarNear(mean, T(1), 0.1));
        expect<RandomSource>(var > std::numeric_limits<T>::epsilon());
    }
}

int main() {
    empty();
    initialize();
    return 0;
}
