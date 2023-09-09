/*
 * Copyright 2023 WeiBo He.
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
#include <fstream>
#include <gperftools/profiler.h>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/TRPMDThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.cuh"
#include "Physica/Core/Parallel/Executor/CudaExecutor.cuh"
#include "Physica/Core/Math/Random/RandomPool.h"

using namespace Physica::Core;
using ScalarType = Scalar<Float>;
using PosScalarType = ScalarType;
using HostForceModel = SilveraGoldman<ScalarType, PosScalarType>;
using DeviceForceModel = device_obj<HostForceModel>;
using RandomPoolType = RandomPool<std::mt19937, 10000>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr unsigned int numMolecular = 108;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

template<class RandomGenerator>
RPMD<ScalarType, PosScalarType> makeSystem(RandomGenerator& gen) {
    using MDCellType = typename RPMD<ScalarType, PosScalarType>::MDCellType;
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    std::uniform_real_distribution dist{};
    for (auto& elem : pos)
        elem = dist(gen);
    typename MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return RPMD<ScalarType, PosScalarType>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}
/**
 * Reference:
 * [1] Miller TF, Manolopoulos DE. 2005. Quantum diffusion in liquid para-hydrogen from ring polymer molecular dynamics. J. Chem. Phys. 122:184503
 */
int main() {
    auto& gen = RandomPoolType::getGen();
    RPMD<ScalarType, PosScalarType> rpmd = makeSystem(gen);
    rpmd.initMomentum(gen);

    HostForceModel hostModel(pair_cutoff);
    DeviceForceModel deviceModel(numMolecular, pair_cutoff);
    const auto f0 = hostModel.template force<SequentialExecutor, true>(rpmd.phaseToCell(0));
    const auto f1 = deviceModel.template force<CudaExecutor, true>(rpmd.phaseToCell(0));
    if (!vectorNear(f0, f1, 1E-4))
        return 1;
    return 0;
}
