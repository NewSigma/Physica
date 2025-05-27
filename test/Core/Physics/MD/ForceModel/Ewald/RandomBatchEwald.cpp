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
#include <iostream>
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/ForceModel/Ewald/RandomBatchEwald.h"

using namespace Physica;
using ScalarType = float64;
using RandomSource = Random<MT19937, 12989825518855205292UL>;
using ForceModel = Q_TIP4P<ScalarType, RandomBatchEwald<ScalarType, RandomSource>>;
using KineticModel = FreeModel<ScalarType, 3, Physica::Dynamic, RPMDIntegrator::Exact>;
using ThermoType = DoubleThermo<KineticModel>;
constexpr size_t numReplica = 32;
constexpr size_t numContract = 8;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(298);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.25;
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(9);
constexpr double massMoleculeInSI = PhyConst<SI>::atomMass(1) * 2 + PhyConst<SI>::atomMass(8);

Vector3D<ScalarType> randomVector() {
    const ScalarType theta(ScalarType::random_uniform<RandomSource>() * M_PI);
    const ScalarType phi(ScalarType::random_uniform<RandomSource>() * M_PI * 2);
    Vector3D<ScalarType> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
    result *= ScalarType(ForceModel::equalR);
    return result;
}

MDCell<ScalarType> makeSystem(unsigned int cellSize) {
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
    auto& gen = RandomSource::getInstance();
    for (size_t i = 0; i < MoleculePerCell; ++i) {
        auto posO = pos.row(i + maxIndexH);
        if (i == 0) {
            posO[0] = latticeFactor * (0.25 + (dist(gen) - 0.5) / 5);
            posO[1] = latticeFactor * (0.25 + (dist(gen) - 0.5) / 5);
            posO[2] = latticeFactor * (0.25 + (dist(gen) - 0.5) / 5);
        }
        else if (i == 1) {
            posO[0] = latticeFactor * (0.75 + (dist(gen) - 0.5) / 5);
            posO[1] = latticeFactor * (0.75 + (dist(gen) - 0.5) / 5);
            posO[2] = latticeFactor * (0.25 + (dist(gen) - 0.5) / 5);
        }
        else if (i == 2) {
            posO[0] = latticeFactor * (0.75 + (dist(gen) - 0.5) / 5);
            posO[1] = latticeFactor * (0.25 + (dist(gen) - 0.5) / 5);
            posO[2] = latticeFactor * (0.75 + (dist(gen) - 0.5) / 5);
        }
        else if (i == 3) {
            posO[0] = latticeFactor * (0.25 + (dist(gen) - 0.5) / 5);
            posO[1] = latticeFactor * (0.75 + (dist(gen) - 0.5) / 5);
            posO[2] = latticeFactor * (0.75 + (dist(gen) - 0.5) / 5);
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
/**
 * Reference:
 * [1] J. Chem. Phys. 131, 024501 (2009); https://doi.org/10.1063/1.3167790
 */
int main() {
    auto cell = makeSystem(2);
    ForceModel::sortPosition(cell);
    ForceModel forceModel(cell, pair_cutoff, RandomBatchEwald<ScalarType, RandomSource>(1000, 200));
    RPMD<ScalarType> rpmd(std::move(cell), numReplica, numContract, temperatureT, timeStep);
    rpmd.initMomentum<KineticModel, RandomSource>();

    constexpr double answer = PhyConst<AU>::angstormToBohr(0.978);
    ScalarType bond = 0;

    ThreadPool::numThreadRequired = 4;
    {
        const ThermoType thermo(temperatureT, thermostatTime);
        KineticModel kineticModel(temperatureT, numReplica);
        rpmd.nvt_step_for<ThermoType, RandomSource, KineticModel, decltype(forceModel), Sequential>(
            PhyConst<AU>::secondToTime(1 * 1E-12),
            thermo,
            kineticModel,
            forceModel);
        for (size_t i = 0; i < 100; ++i) {
            const PeriodicCell<ScalarType, 3> cell = rpmd.makeAverageCell();
            ScalarType temp = 0;
            const size_t numO = rpmd.getNumParticle() / 3;
            const size_t numH = numO * 2;
            for (size_t j = 0; j < numO; ++j) {
                toNextMean(temp, 2 * j, cell.minDistVector(numH + j, 2 * j).norm());
                toNextMean(temp, 2 * j + 1, cell.minDistVector(numH + j, 2 * j + 1).norm());
            }
            toNextMean(bond, i, temp);
            rpmd.nvt_step<ThermoType, RandomSource, KineticModel, decltype(forceModel), Sequential>(thermo, kineticModel, forceModel);
        }
    }
    ThreadPool::getInstance().shouldExit();
    if (!scalarNear(bond, ScalarType(answer), 2E-2))
        return 1;
    return 0;
}
