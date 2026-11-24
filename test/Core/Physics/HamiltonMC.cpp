/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
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
        hmc.setTimeStep(T(1));
        hmc.step<RandomSource>(EmptyForceModel<T, 1>{});
    }

    void initialize() {
        // Test that temperature is not fixed
        using KineticModel = OpenModel<T, 1, 1>;
        constexpr size_t numSample = 16384; // standard error is 1 / sqrt(numSample); P(fail) < 1E-9
        auto hmc = HamiltonMC<T>({1, 1});
        auto& root = hmc.getRoot();
        T mean = 0;
        T var = 0;
        for (size_t i = 0; i < numSample; ++i) {
            root.template initMomentum<KineticModel, RandomSource>();
            var.toNextVariance(mean, i, root.template calcTemperature<KineticModel>());
        }
        expect<RandomSource>(scalarNear(mean, T(1), 0.05));
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
