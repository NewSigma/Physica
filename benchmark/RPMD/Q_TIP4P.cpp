/*
 * Copyright 2022-2024 Weibo He.
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
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include <benchmark/benchmark.h>
#include <gperftools/profiler.h>
#include <iostream>
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Utils/BenchmarkHelper.h"

using namespace Physica;
using ScalarType = float64;
using KineticModel = FreeModel<ScalarType, 3, Physica::Dynamic, RPMDIntegrator::Exact>;
using ForceModel = Q_TIP4P<ScalarType, Ewald<ScalarType>>;
using ThermoType = DoubleThermo<KineticModel>;
using RandomType = Random<MT19937>;
constexpr size_t numReplica = 32;
constexpr size_t numContract = 8;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(298);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.1;
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(9);
constexpr double massMoleculeInSI = PhyConst<SI>::atomMass(1) * 2 + PhyConst<SI>::atomMass(8);

static Vector3D<ScalarType> randomVector() {
    const ScalarType theta(ScalarType::random_uniform<RandomType>() * M_PI);
    const ScalarType phi(ScalarType::random_uniform<RandomType>() * M_PI * 2);
    Vector3D<ScalarType> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
    result *= ScalarType(ForceModel::equalR);
    return result;
}

static MDCell<ScalarType> makeSystem(unsigned int cellSize) {
    using CrystalCellType = CrystalCell<ScalarType>;
    constexpr size_t MoleculePerCell = 4;
    constexpr size_t maxIndexH = MoleculePerCell * 2;
    constexpr size_t maxIndexO = MoleculePerCell * 3;
    constexpr size_t numAtom = MoleculePerCell * 3;

    ScalarType cellVolume = ((MoleculePerCell * massMoleculeInSI * 1000 / 0.997) * 1E-6) / (PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius);
    const ScalarType latticeFactor(cbrt(cellVolume));
    CrystalCellType::LatticeMatrix lattice = CrystalCellType::LatticeMatrix::unitMatrix(3);
    lattice *= latticeFactor;

    CrystalCellType::PositionMatrix pos(numAtom, 3);
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
        posH1 = posO + randomVector();
        posH2 = posO + randomVector();
    }

    CrystalCellType::AtomicArray atomicNumbers(numAtom);
    for (size_t i = 0; i < maxIndexH; ++i)
        atomicNumbers[i] = 1;
    for (size_t i = maxIndexH; i < maxIndexO; ++i)
        atomicNumbers[i] = 8;

    CrystalCellType cell({std::move(lattice), std::move(pos), CrystalCellType::Type::Cartesian}, std::move(atomicNumbers));
    cell.toSuperCell(cellSize, cellSize, cellSize);
    return MDCell<ScalarType>(std::move(cell));
}

static void bench(benchmark::State& state) {
    auto cell = makeSystem(2);
    ForceModel::sortPosition(cell);
    ForceModel forceModel(cell, pair_cutoff, {});
    KineticModel kineticModel(temperatureT, numReplica);
    const ThermoType thermo(temperatureT, thermostatTime);
    RPMD<ScalarType> rpmd(std::move(cell), numReplica, numContract, temperatureT, timeStep);
    rpmd.initMomentum<KineticModel, RandomType>();
    for (auto _ : state)
        rpmd.nve_step<KineticModel, ForceModel, SeqExecutor>(kineticModel, forceModel);
}

BENCHMARK(bench)->Name("Q_TIP4P")->Unit(benchmark::kMillisecond);
