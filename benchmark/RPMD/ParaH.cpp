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
#include <gperftools/profiler.h>
#include <iostream>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/ForceModel/SilveraGoldman.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Utils/BenchmarkHelper.h"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Float>;
using KineticModel = FreeModel<ScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
using ForceModel = SilveraGoldman<ScalarType>;
using RandomPoolType = RandomPool<std::mt19937, 3438603950906262893>;
constexpr size_t numReplica = 24;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(25);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
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
    ThreadPool::numThreadRequired = 4;
    ThreadPool& pool = ThreadPool::getInstance();

    auto& gen = RandomPoolType::getGen();
    KineticModel kineticModel(temperatureT, numReplica);
    ForceModel forceModel(pair_cutoff);
    {
        auto rpmd = makeSystem(108, gen);
        rpmd.initMomentum(gen);
        rpmd.updateForce<ForceModel, ThreadExecutor>(forceModel);
        auto timeuse = Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-13),
                kineticModel,
                forceModel);
        }, 8, 20);
        std::cout << "108 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        auto rpmd = makeSystem(256, gen);
        rpmd.initMomentum(gen);
        rpmd.updateForce<ForceModel, ThreadExecutor>(forceModel);
        auto timeuse = Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-13),
                kineticModel,
                forceModel);
        }, 8, 20);
        std::cout << "256 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        auto rpmd = makeSystem(500, gen);
        rpmd.initMomentum(gen);
        rpmd.updateForce<ForceModel, ThreadExecutor>(forceModel);
        auto timeuse = Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(1 * 1E-13),
                kineticModel,
                forceModel);
        }, 8, 20);
        std::cout << "500 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        auto rpmd = makeSystem(864, gen);
        rpmd.initMomentum(gen);
        rpmd.updateForce<ForceModel, ThreadExecutor>(forceModel);
        auto timeuse = Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(5 * 1E-14),
                kineticModel,
                forceModel);
        }, 8, 20);
        std::cout << "864 atom time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    pool.shouldExit();
    return 0;
}
