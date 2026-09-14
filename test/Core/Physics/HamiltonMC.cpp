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

    class InfiniteWallModel {
        size_t numCall = 0;
    public:
        template<ExecutePolicy P>
        [[nodiscard]] T potentialV(const MDCell<T, 1>&) {
            return (numCall++ % 2 == 0) ? T(1E300) : T(0);
        }
    };

    void infinity_accept() {
        // Test that we do not trigger a exp overflow if energy decrease is large
        auto hmc = HamiltonMC<T>({1, 1});
        auto& root = hmc.getRoot();
        root.getPhaseMatrix().col(0).tail(hmc.getDOF()) = T(1);

        InfiniteWallModel wall{};
        const bool accept = hmc.step_radial<RandomSource>(wall, T(1));
        expect<RandomSource>(accept);
        expect<RandomSource>(hmc.getSample().isFinite());
    }
}

int main() {
    empty();
    initialize();
    infinity_accept();
    return 0;
}
