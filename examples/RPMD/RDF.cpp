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
#include <fstream>
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/ElectronicStructure/CrystalCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Utils/Random.h"
#include "Physica/Core/Physics/MD/RDF.h"
#include "Physica/Core/IO/Poscar.h"

using namespace Physica::Core;
using namespace Physica::Core::Parallel;
using ScalarType = Scalar<Double, false>;
using PosScalarType = Scalar<Double, false>;
using ThermostatType = DoubleThermo<ScalarType, PosScalarType>;
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

    CrystalCell cell({std::move(lattice), std::move(pos), CrystalCell::Type::Cartesian}, std::move(atomicNumbers));
    cell.unitToSuper(cellSize, cellSize, cellSize);
    return MDCell<ScalarType, PosScalarType>(std::move(cell));
}
/**
 * Reference:
 * [1] S. Habershon, T. E. Markland, and D. E. Manolopoulosa, J. Chem. Phys. 131, 024501(2009)
 */
int main() {
    std::mt19937 gen{};

    auto cell = makeSystem(3, gen);
    ForceModel::sortPosition(cell);
    const ThermostatType thermo(temperatureT, thermostatTime);
    RPMD<ScalarType, PosScalarType> rpmd(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
    rpmd.initMomentum(gen);
    ForceModel model(rpmd.phaseToCell(0), pair_cutoff);

    RDF<PosScalarType> rdf;
    {
        Physica::Utils::Array<bool> isFromParticle(rpmd.getNumParticle());
        Physica::Utils::Array<bool> isToParticle(rpmd.getNumParticle());
        for (size_t i = 0; i < isFromParticle.getLength(); ++i) {
            const bool isHydrogen = i < isFromParticle.getLength() * 2 / 3;
            isFromParticle[i] = isHydrogen;
            isToParticle[i] = isHydrogen;
        }
        rdf = RDF<PosScalarType>(std::move(isFromParticle), std::move(isToParticle), PhyConst<AU>::angstormToBohr(0.01), 700);
    }

    ThreadPool::initThreadPool(4);
    ThreadPool& pool = ThreadPool::getInstance();
    {
        rpmd.nvt_step_for<ThermostatType, decltype(gen), decltype(model), ThreadExecutor>(PhyConst<AU>::secondToTime(2 * 1E-12), thermo, gen, model);
        for (size_t i = 0; i < 1000; ++i) {
            rpmd.nvt_step<ThermostatType, decltype(gen), decltype(model), ThreadExecutor>(thermo, gen, model);
            rpmd.normalizeCentroid();
            for (size_t j = 0; j < numReplica; ++j)
                rdf.sample(rpmd.phaseToCell(j));
        }
    }
    pool.shouldExit();
    ThreadPool::deInitThreadPool();
    {
        std::ofstream fout("dists");
        fout << rdf.makeDists();
    }
    {
        std::ofstream fout("rdf");
        fout << rdf.makeRDF(rpmd.getVolume());
    }
    {
        Poscar poscar({rpmd.getLattice(), rpmd.getRingPolymer().makeCentroidPos(), Poscar::Type::Cartesian}, {rpmd.getNumParticle() * 2 / 3, rpmd.getNumParticle() / 3});
        std::ofstream fout("H2O.vasp");
        fout << poscar;
        fout.close();
    }
    return 0;
}
