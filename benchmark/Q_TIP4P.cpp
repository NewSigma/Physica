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
#include <gperftools/profiler.h>
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Utils/Cycler.h"

using namespace Physica::Core;
using namespace Physica::Core::Parallel;
using namespace Physica::Utils;
using ScalarType = Scalar<Double, false>;
using PosScalarType = Scalar<Double, false>;
using ForceModel = Q_TIP4P<ScalarType, PosScalarType, 1>;
constexpr size_t numReplica = 32;
constexpr double temperatureT = PhyConst<AU>::kToTemperature(298);
constexpr double thermostatTime = PhyConst<AU>::secondToTime(100 * 1E-15);
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.1;
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

    typename CrystalCell::AtomicArray atomicNumbers(numAtom);
    for (size_t i = 0; i < maxIndexH; ++i)
        atomicNumbers[i] = 1;
    for (size_t i = maxIndexH; i < maxIndexO; ++i)
        atomicNumbers[i] = 8;

    CrystalCell cell(std::move(lattice), std::move(pos), std::move(atomicNumbers), CrystalCell::Type::Cartesian);
    cell.unitToSuper(cellSize, cellSize, cellSize);
    return MDCell<ScalarType, PosScalarType>(std::move(cell));
}

int main() {
    Cycler::init();
    std::mt19937 gen{};

    RPMD<ScalarType, PosScalarType> rpmd(makeSystem(1, gen), numReplica, numReplica, temperatureT, thermostatTime, timeStep);
    rpmd.initMomentum(gen);
    ForceModel model(rpmd.phaseToCell(0), pair_cutoff);
    {
        const auto from = Cycler::tic();
        rpmd.nvt_step_for<decltype(gen), decltype(model), SequentialExecutor>(PhyConst<AU>::secondToTime(1 * 1E-13), gen, model);
        const auto to = Cycler::toc();
        std::cout << "1 Threads time use: " << Cycler::toSeconds(to - from) << '\n';
    }
    {
        ThreadPool::initThreadPool(2);
        ThreadPool& pool = ThreadPool::getInstance();
        const auto from = Cycler::tic();
        rpmd.nvt_step_for<decltype(gen), decltype(model), ThreadExecutor>(PhyConst<AU>::secondToTime(1 * 1E-13), gen, model);
        const auto to = Cycler::toc();
        std::cout << "2 Threads time use: " << Cycler::toSeconds(to - from) << '\n';
        pool.shouldExit();
        ThreadPool::deInitThreadPool();
    }
    {
        ThreadPool::initThreadPool(4);
        ThreadPool& pool = ThreadPool::getInstance();
        const auto from = Cycler::tic();
        rpmd.nvt_step_for<decltype(gen), decltype(model), ThreadExecutor>(PhyConst<AU>::secondToTime(1 * 1E-13), gen, model);
        const auto to = Cycler::toc();
        std::cout << "4 Threads time use: " << Cycler::toSeconds(to - from) << '\n';
        pool.shouldExit();
        ThreadPool::deInitThreadPool();
    }
    return 0;
}
