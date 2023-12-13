/*
 * Copyright 2022 WeiBo He.
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
#include <gperftools/profiler.h>
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Utils/BenchmarkHelper.h"

using namespace Physica::Core;
using namespace Physica::Utils;
using ScalarType = Scalar<Double>;
using ThermostatType = DoubleThermo<ScalarType>;
using KineticModel = FreeModel<ScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
using ForceModel = Q_TIP4P<ScalarType, Ewald<ScalarType>, false>;
using RandomPoolType = RandomPool<std::mt19937>;
constexpr size_t numReplica = 32;
constexpr size_t numContract = 8;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(298);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.1;
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(9);
constexpr double massMoleculeInSI = PhyConst<SI>::atomMass(1) * 2 + PhyConst<SI>::atomMass(8);

template<class RandomGenerator>
Vector<ScalarType, 3> randomVector(RandomGenerator& gen) {
    std::uniform_real_distribution dist{};
    const ScalarType theta(dist(gen) * M_PI);
    const ScalarType phi(dist(gen) * M_PI * 2);
    Vector<ScalarType, 3> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
    result *= ScalarType(ForceModel::equalR);
    return result;
}

template<class RandomGenerator>
MDCell<ScalarType> makeSystem(unsigned int cellSize, RandomGenerator& gen) {
    using CrystalCellType = CrystalCell<ScalarType>;
    constexpr size_t MoleculePerCell = 4;
    constexpr size_t maxIndexH = MoleculePerCell * 2;
    constexpr size_t maxIndexO = MoleculePerCell * 3;
    constexpr size_t numAtom = MoleculePerCell * 3;

    ScalarType cellVolume = ((MoleculePerCell * massMoleculeInSI * 1000 / 0.997) * 1E-6) / (PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius);
    const ScalarType latticeFactor(cbrt(cellVolume));
    typename CrystalCellType::LatticeMatrix lattice = CrystalCellType::LatticeMatrix::unitMatrix(3);
    lattice *= latticeFactor;

    typename CrystalCellType::PositionMatrix pos(numAtom, 3);
    std::uniform_real_distribution dist{};
    for (size_t i = 0; i < MoleculePerCell; ++i) {
        auto posO = pos.row(i + maxIndexH);
        if (i == 0) {
            posO[0] = 0;
            posO[1] = 0;
            posO[2] = 0;
        }
        if (i == 1) {
            posO[0] = latticeFactor * 0.5;
            posO[1] = latticeFactor * 0.5;
            posO[2] = 0;
        }
        if (i == 2) {
            posO[0] = latticeFactor * 0.5;
            posO[1] = 0;
            posO[2] = latticeFactor * 0.5;
        }
        if (i == 3) {
            posO[0] = 0;
            posO[1] = latticeFactor * 0.5;
            posO[2] = latticeFactor * 0.5;
        }
        auto posH1 = pos.row(2 * i);
        auto posH2 = pos.row(2 * i + 1);
        posH1 = posO + randomVector(gen);
        posH2 = posO + randomVector(gen);
    }

    typename CrystalCellType::AtomicArray atomicNumbers(numAtom);
    for (size_t i = 0; i < maxIndexH; ++i)
        atomicNumbers[i] = 1;
    for (size_t i = maxIndexH; i < maxIndexO; ++i)
        atomicNumbers[i] = 8;

    CrystalCellType cell({std::move(lattice), std::move(pos), CrystalCellType::Type::Cartesian}, std::move(atomicNumbers));
    cell.toSuperCell(cellSize, cellSize, cellSize);
    return MDCell<ScalarType>(std::move(cell));
}

int main() {
    auto& gen = RandomPoolType::getGen();
    auto cell = makeSystem(2, gen);
    ForceModel::sortPosition(cell);
    ForceModel forceModel(cell, pair_cutoff, {});
    KineticModel kineticModel(temperatureT, numReplica);
    const ThermostatType thermo(temperatureT, thermostatTime);
    RPMD<ScalarType> rpmd(std::move(cell), numReplica, numContract, temperatureT, timeStep);
    rpmd.initMomentum(gen);
    {
        auto timeuse = Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, SequentialExecutor>(
                PhyConst<AU>::secondToTime(1 * 1E-14),
                kineticModel,
                forceModel);
        }, 8, 20);
        std::cout << "1 thread time use: " << timeuse.first << '(' << timeuse.second << ")\n";
    }
    {
        ThreadPool::numThreadRequired = 2;
        auto timeuse = Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-14),
                kineticModel,
                forceModel);
        }, 8, 20);
        std::cout << "2 thread time use: " << timeuse.first << '(' << timeuse.second << ")\n";
        ThreadPool::getInstance().shouldExit();
    }
    {
        ThreadPool::numThreadRequired = 4;
        ThreadPool& pool = ThreadPool::getInstance();
        pool.restart();
        auto timeuse = Benchmark::run([&]() {
            rpmd.nve_step_for<KineticModel, ForceModel, ThreadExecutor>(
                PhyConst<AU>::secondToTime(2 * 1E-14),
                kineticModel,
                forceModel);
        }, 8, 20);
        std::cout << "4 thread time use: " << timeuse.first << '(' << timeuse.second << ")\n";
        pool.shouldExit();
    }
    return 0;
}
