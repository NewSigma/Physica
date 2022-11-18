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
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Utils/Random.h"

using namespace Physica::Core;
using namespace Physica::Core::Parallel;
using ScalarType = Scalar<Double, false>;
using PosScalarType = Scalar<Double, false>;
using ForceModel = Q_TIP4P<ScalarType, PosScalarType>;
constexpr size_t numReplica = 32;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(298);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.25;
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(9);
constexpr double massMoleculeInSI = PhyConst<SI>::atomMass(1) * 2 + PhyConst<SI>::atomMass(8);

template<class RandomGenerator>
Vector<PosScalarType, 3> randomVector(RandomGenerator& gen) {
    std::uniform_real_distribution dist{};
    const PosScalarType theta(dist(gen) * M_PI);
    const PosScalarType phi(dist(gen) * M_PI * 2);
    Vector<PosScalarType, 3> result{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
    result *= PosScalarType(ForceModel::equalR);
    return result;
}

template<class RandomGenerator>
MDCell<ScalarType, PosScalarType> makeSystem(unsigned int cellSize, RandomGenerator& gen) {
    constexpr size_t MoleculePerCell = 4;
    constexpr size_t maxIndexH = MoleculePerCell * 2;
    constexpr size_t maxIndexO = MoleculePerCell * 3;
    constexpr size_t numAtom = MoleculePerCell * 3;

    PosScalarType cellVolume = ((MoleculePerCell * massMoleculeInSI * 1000 / 0.997) * 1E-6) / (PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius * PhyConst<SI>::bohrRadius);
    const PosScalarType latticeFactor(cbrt(cellVolume));
    typename CrystalCell::LatticeMatrix lattice = CrystalCell::LatticeMatrix::unitMatrix(3);
    lattice *= latticeFactor;

    typename CrystalCell::PositionMatrix pos(numAtom, 3);
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

    typename CrystalCell::AtomicArray atomicNumbers(numAtom);
    for (size_t i = 0; i < maxIndexH; ++i)
        atomicNumbers[i] = 1;
    for (size_t i = maxIndexH; i < maxIndexO; ++i)
        atomicNumbers[i] = 8;

    CrystalCell cell(std::move(lattice), std::move(pos), std::move(atomicNumbers), CrystalCell::Type::Cartesian);
    cell.unitToSuper(cellSize, cellSize, cellSize);
    return MDCell<ScalarType, PosScalarType>(std::move(cell));
}

void testSort() {
    using CellType = MDCell<ScalarType, PosScalarType>;
    using PositionMatrix = typename CellType::PositionMatrix;
    std::mt19937 gen{};
    CellType cell = makeSystem(3, gen);
    const CellType origin_cell = cell;
    auto order = ForceModel::sortPosition(cell);
    PositionMatrix result = order.transpose() * cell.getPos();
    if (!matrixNear(result, origin_cell.getPos(), 1E-15))
        exit(EXIT_FAILURE);
}

void testForce() {
    std::mt19937 gen(9806048078107704755UL);
    auto cell = makeSystem(2, gen);
    ForceModel::sortPosition(cell);
    ForceModel model(cell, pair_cutoff);
    const auto& pos = cell.getPos();
    const Vector<ScalarType> force1 = model.template force<SequentialExecutor>(cell);
    Vector<ScalarType> force2(force1.getLength());
    /* Force from finite differential */ {
        for (size_t i = 0; i < force2.getLength(); ++i) {
            force2[i] = -Differential<ScalarType>::ridders([i, &cell, &model](ScalarType x) -> ScalarType {
                typename MDCell<ScalarType, PosScalarType>::PositionMatrix temp = cell.getPos();
                *(temp.begin() + i) = PosScalarType(x);
                MDCell<ScalarType, PosScalarType> new_cell(cell.getLattice(), std::move(temp), cell.getMassVec());
                return model.potentialEnergy(new_cell);
            }, pos.flatten().calc(i), 0.3);
        }
    }
    if (!vectorNear(force1, force2, 1E-3))
        exit(EXIT_FAILURE);
}

void testMD() {
    std::mt19937 gen(12989825518855205292UL);

    auto cell = makeSystem(2, gen);
    ForceModel::sortPosition(cell);
    ForceModel model(cell, pair_cutoff);
    RPMD<ScalarType, PosScalarType> rpmd(std::move(cell), numReplica, numReplica, temperatureT, thermostatTime, timeStep);
    rpmd.initMomentum(gen);

    constexpr double answer = PhyConst<AU>::angstormToBohr(0.978);
    ScalarType bond = 0;

    ThreadPool::initThreadPool(4);
    ThreadPool& pool = ThreadPool::getInstance();
    {
        rpmd.nvt_step_for<decltype(gen), decltype(model), ThreadExecutor>(PhyConst<AU>::secondToTime(1 * 1E-12), gen, model);
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
            rpmd.nvt_step<decltype(gen), decltype(model), ThreadExecutor>(gen, model);
        }
    }
    pool.shouldExit();
    ThreadPool::deInitThreadPool();
    if (!scalarNear(bond, ScalarType(answer), 2E-2))
        exit(EXIT_FAILURE);
}
/**
 * Reference:
 * [1] S. Habershon, T. E. Markland, and D. E. Manolopoulosa, J. Chem. Phys. 131, 024501(2009)
 */
int main() {
    //testSort();
    //testForce();
    testMD();
    return 0;
}
