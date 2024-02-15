/*
 * Copyright 2023-2024 WeiBo He.
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
#include "Physica/Core/Physics/MD/ForceModel/CPUGPUModel.cuh"
#include "Physica/Core/Math/Random/RandomPool.h"
#include "Physica/Utils/BenchmarkHelper.h"

using namespace Physica::Core;
using ScalarType = Scalar<Float>;
using MDType = RPMD<ScalarType, 3, Physica::Utils::Dynamic, Physica::Utils::PageLockedAllocator<ScalarType>>;
using MDCellType = typename MDType::MDCellType;
using KineticModel = FreeModel<ScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
using RandomPoolType = RandomPool<std::mt19937, 10000>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr double pair_cutoff = 15;
constexpr double molarVolume = 31.7;
constexpr double mass = PhyConst<AU>::atomMass(1) * 2;

template<class RandomGenerator>
MDType makeSystem(size_t numMolecular, RandomGenerator& gen) {
    typename MDCellType::LatticeMatrix lattice = MDCellType::LatticeMatrix::unitMatrix(3);
    typename MDCellType::PositionMatrix pos(numMolecular, 3);
    std::uniform_real_distribution dist{};
    for (auto& elem : pos)
        elem = dist(gen);
    typename MDCellType::MassVector massVec(numMolecular, mass);
    MDCellType cell(std::move(lattice), std::move(pos), std::move(massVec));

    const double factor = (std::cbrt(numMolecular * molarVolume / PhyConst<SI>::avogadroNa) / 100) / PhyConst<SI>::bohrRadius;
    cell.scale(factor);

    return MDType(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
}
/**
 * Reference:
 * [1] Miller TF, Manolopoulos DE. 2005. Quantum diffusion in liquid para-hydrogen from ring polymer molecular dynamics. J. Chem. Phys. 122:184503
 */
int main() {
    ThreadPool::numThreadRequired = 8;
    auto& gen = RandomPoolType::getInstance().getGen();
    KineticModel kineticModel(temperatureT, numReplica);
    {
        using HostModel = SilveraGoldman<ScalarType, true, true>;
        using ForceModel = CPUGPUModel<HostModel>;
        constexpr size_t numMolecular = 108;
        MDType rpmd = makeSystem(numMolecular, gen);
        rpmd.initMomentum(gen);
        ForceModel forceModel(4, HostModel(pair_cutoff), numMolecular);
        auto timeuse = Physica::Utils::Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, AutoExecutor>(PhyConst<AU>::secondToTime(2 * 1E-13), kineticModel, forceModel);
        }, 8, 20);
        std::cout << "108 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        using HostModel = SilveraGoldman<ScalarType, true, true>;
        using ForceModel = CPUGPUModel<HostModel>;
        constexpr size_t numMolecular = 256;
        MDType rpmd = makeSystem(numMolecular, gen);
        rpmd.initMomentum(gen);
        ForceModel forceModel(4, HostModel(pair_cutoff), numMolecular);
        auto timeuse = Physica::Utils::Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, AutoExecutor>(PhyConst<AU>::secondToTime(2 * 1E-13), kineticModel, forceModel);
        }, 8, 20);
        std::cout << "256 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        using HostModel = SilveraGoldman<ScalarType, true, true>;
        using ForceModel = CPUGPUModel<HostModel>;
        constexpr size_t numMolecular = 500;
        MDType rpmd = makeSystem(numMolecular, gen);
        rpmd.initMomentum(gen);
        ForceModel forceModel(4, HostModel(pair_cutoff), numMolecular);
        auto timeuse = Physica::Utils::Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, AutoExecutor>(PhyConst<AU>::secondToTime(1 * 1E-13), kineticModel, forceModel);
        }, 8, 20);
        std::cout << "500 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        using HostModel = SilveraGoldman<ScalarType, true, false>;
        using ForceModel = CPUGPUModel<HostModel>;
        constexpr size_t numMolecular = 864;
        MDType rpmd = makeSystem(numMolecular, gen);
        rpmd.initMomentum(gen);
        ForceModel forceModel(4, HostModel(pair_cutoff), numMolecular);
        auto timeuse = Physica::Utils::Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, AutoExecutor>(PhyConst<AU>::secondToTime(5 * 1E-14), kineticModel, forceModel);
        }, 8, 20);
        std::cout << "864 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    return 0;
}
