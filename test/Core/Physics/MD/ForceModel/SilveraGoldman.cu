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
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.cuh"
#include "Physica/Core/Parallel/Executor/CUDAExecutor.cuh"
#include "Physica/Core/Math/Random/Random.h"

using namespace Physica;
using ScalarType = float32;
using HostForceModel = SilveraGoldman<ScalarType, true, true>;
using DeviceForceModel = device_obj<HostForceModel>;
using RandomSource = Random<MT19937, 15522090741289029828UL>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

RPMD<ScalarType> makeSystem(size_t numMolecular) {
    using MDCellType = RPMD<ScalarType>::MDCellType;
    MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    auto pos = MDCellType::PositionMatrix::random_uniform<RandomSource>(numMolecular, 3);
    MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return RPMD<ScalarType>(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}
/**
 * Reference:
 * [1] J. Chem. Phys. 122, 184503 (2005); https://doi.org/10.1063/1.1893956
 */
int main() {
    HostForceModel hostModel(pair_cutoff);
    for (size_t numMolecular : {108, 160}) {
        DeviceForceModel deviceModel(numMolecular, pair_cutoff);
        RPMD<ScalarType> rpmd = makeSystem(numMolecular);
        const auto f0 = hostModel.template force<SeqExecutor>(rpmd.phaseToCell(0));
        const auto f1 = deviceModel.template force<CUDAExecutor>(rpmd.phaseToCell(0));
        if (!vectorNear(f0, f1, 1E-3))
            return 1;
    }
    {
        constexpr unsigned int numMolecular = 108;
        DeviceForceModel deviceModel(numMolecular, pair_cutoff);
        RPMD<ScalarType> rpmd = makeSystem(numMolecular);
        const auto f0 = hostModel.template force<SeqExecutor>(rpmd.phaseToCell(0));
        const auto f1 = deviceModel.template force<CUDAExecutor>(rpmd.phaseToCell(0));
        if (!vectorNear(f0, f1, 1E-4))
            return 1;
    }
    return 0;
}
