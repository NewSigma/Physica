/*
 * Copyright 2023-2024 Weibo He.
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
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/TRPMDThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.cuh"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Physica/Core/Math/Random/RandomPool.h"

using namespace Physica::Core;
using ScalarType = float32;
using HostForceModel = SilveraGoldman<ScalarType, true, true>;
using DeviceForceModel = device_obj<HostForceModel>;
using RandomPoolType = RandomPool<std::mt19937, 10000>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

template<class RandomGenerator>
RPMD<ScalarType> makeSystem(size_t numMolecular, RandomGenerator& gen) {
    using MDCellType = typename RPMD<ScalarType>::MDCellType;
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    std::uniform_real_distribution dist{};
    for (auto& elem : pos)
        elem = dist(gen);
    typename MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return RPMD<ScalarType>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}
/**
 * Reference:
 * [1] Miller TF, Manolopoulos DE. 2005. Quantum diffusion in liquid para-hydrogen from ring polymer molecular dynamics. J. Chem. Phys. 122:184503
 */
int main() {
    auto& gen = RandomPoolType::getInstance().getGen();
    HostForceModel hostModel(pair_cutoff);
    for (size_t numMolecular : {108, 160}) {
        DeviceForceModel deviceModel(numMolecular, pair_cutoff);
        RPMD<ScalarType> rpmd = makeSystem(numMolecular, gen);
        const auto f0 = hostModel.template force<SequentialExecutor>(rpmd.phaseToCell(0));
        const auto f1 = deviceModel.template force<CUDAExecutor>(rpmd.phaseToCell(0));
        if (!vectorNear(f0, f1, 1E-3))
            return 1;
    }
    {
        constexpr unsigned int numMolecular = 108;
        DeviceForceModel deviceModel(numMolecular, pair_cutoff);
        RPMD<ScalarType> rpmd = makeSystem(numMolecular, gen);
        const auto f0 = hostModel.template force<SequentialExecutor>(rpmd.phaseToCell(0));
        const auto f1 = deviceModel.template force<CUDAExecutor>(rpmd.phaseToCell(0));
        if (!vectorNear(f0, f1, 1E-4))
            return 1;
    }
    return 0;
}
