/*
 * Copyright 2022-2023 WeiBo He.
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
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double>;
using ThermostatType = DoubleThermo<ScalarType>;
using ForceModel = Q_TIP4P<ScalarType, Ewald<ScalarType>>;
using KineticModel = FreeModel<ScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
using RandomPoolType = RandomPool<std::mt19937, 12989825518855205292UL>;
constexpr size_t numReplica = 32;
constexpr size_t numContract = 8;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(298);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.25;
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

void testSort() {
    using CellType = MDCell<ScalarType>;
    using PositionMatrix = typename CellType::PositionMatrix;
    std::mt19937 gen{};
    CellType cell = makeSystem(3, gen);
    const CellType origin_cell = cell;
    auto order = ForceModel::sortPosition(cell);
    PositionMatrix result = order.transpose() * cell.getPos();
    if (!matrixNear(result, origin_cell.getPos(), 1E-15))
        exit(EXIT_FAILURE);
}

void testMD() {
    auto& gen = RandomPoolType::getGen();
    auto cell = makeSystem(2, gen);
    ForceModel::sortPosition(cell);
    ForceModel forceModel(cell, pair_cutoff, Ewald<ScalarType>{});
    RPMD<ScalarType> rpmd(std::move(cell), numReplica, numContract, temperatureT, timeStep);
    rpmd.initMomentum(gen);

    constexpr double answer = PhyConst<AU>::angstormToBohr(0.978);
    ScalarType bond = 0;

    ThreadPool::numThreadRequired = 4;
    {
        const ThermostatType thermo(temperatureT, thermostatTime);
        KineticModel kineticModel(temperatureT, numReplica);
        rpmd.nvt_step_for<ThermostatType, RandomPoolType, KineticModel, decltype(forceModel), ThreadExecutor>(
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
            rpmd.nvt_step<ThermostatType, RandomPoolType, KineticModel, decltype(forceModel), ThreadExecutor>(thermo, kineticModel, forceModel);
        }
    }
    ThreadPool::getInstance().shouldExit();
    if (!scalarNear(bond, ScalarType(answer), 2E-2))
        exit(EXIT_FAILURE);
}
/**
 * Reference:
 * [1] S. Habershon, T. E. Markland, and D. E. Manolopoulosa, J. Chem. Phys. 131, 024501(2009)
 */
int main() {
    testSort();
    testMD();
    return 0;
}
