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
#include "Physica/Core/Physics/ManyBody/FreqDQMC.h"
#include "Physica/Core/Physics/MC/HamiltonMC.h"
#include "Test.h"

using namespace Physica;
using T = float64;
using Tc = cfloat64;
using RandomSource = Random<PCG64DXSM>;
constexpr T HoppingT = 0;
constexpr T RepelU = 8;
constexpr T Beta = 8;
constexpr int Dim = 1;
constexpr int NumSite = 1;
constexpr int FreqDensity = 0;
constexpr T StepSize = 1E-3;
constexpr T Duration = 10;

int main() {
    const SquareLattice<Dim> lattice({NumSite}, 1);
    const HubbardParams<T> params(HoppingT, RepelU, lattice, Beta, RepelU * 0.5, 1);
    auto dqmc = FreqDQMC<Tc>(params, FreqDensity);

    auto mass = VectorND<T>(dqmc.getAuxField().getSize() * 2, 1);
    auto hmc = HamiltonMC<T>(std::move(mass));
    auto& engine = hmc.getRoot();
    engine.setTimeStep(StepSize);

    using Kinetic = OpenModel<T, 1, 1>;
    auto kinetic = Kinetic(1, 1);
    engine.template initMomentum<Kinetic, RandomSource>();

    const T prevE = engine.calcClassicalInternalEnergy(dqmc);
    engine.nve_step_for(Duration, kinetic, dqmc);
    const T curE = engine.calcClassicalInternalEnergy(dqmc);
    expect(scalarNear(prevE, curE, 1E-4)); // Energe conserves
    return 0;
}
