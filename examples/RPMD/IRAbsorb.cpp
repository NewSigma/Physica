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
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/IRAbsorb.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/IO/Poscar.h"

using namespace Physica::Core;
using ScalarType = Scalar<Double>;
using PosScalarType = Scalar<Double>;
using ThermostatType = DoubleThermo<ScalarType, PosScalarType>;
using ForceModel = Q_TIP4P<ScalarType, PosScalarType>;
using KineticModel = FreeModel<ScalarType, PosScalarType, 3, Dynamic, RPMDIntegrator::Exact>;
using RandomPoolType = RandomPool<std::mt19937, 16172479642808547241UL>;
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
    typename CrystalCell<PosScalarType>::LatticeMatrix lattice = CrystalCell<PosScalarType>::LatticeMatrix::unitMatrix(3);
    lattice *= latticeFactor;

    typename CrystalCell<PosScalarType>::PositionMatrix pos(numAtom, 3);
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

    typename CrystalCell<PosScalarType>::AtomicArray atomicNumbers(numAtom);
    for (size_t i = 0; i < maxIndexH; ++i)
        atomicNumbers[i] = 1;
    for (size_t i = maxIndexH; i < maxIndexO; ++i)
        atomicNumbers[i] = 8;

    CrystalCell<PosScalarType> cell({std::move(lattice), std::move(pos), CrystalCell<PosScalarType>::Type::Cartesian}, std::move(atomicNumbers));
    cell.toSuperCell(cellSize, cellSize, cellSize);
    return MDCell<ScalarType, PosScalarType>(std::move(cell));
}

void testMD() {
    auto& gen = RandomPoolType::getGen();
    auto cell = makeSystem(3, gen);
    ForceModel::sortPosition(cell);
    ThermostatType thermo(temperatureT, thermostatTime);
    RPMD<ScalarType, PosScalarType> rpmd(std::move(cell), numReplica, 8, temperatureT, timeStep);
    rpmd.initMomentum(gen);
    ForceModel forceModel(rpmd.phaseToCell(0), pair_cutoff);

    Vector<ScalarType> corr(40000);

    ThreadPool::numThreadRequired = 16;
    {
        KineticModel kineticModel(temperatureT, numReplica);
        for (size_t path = 0; path < 16; ++path) {
            rpmd.nvt_step_for<ThermostatType, RandomPoolType, KineticModel, ForceModel, ThreadExecutor>(PhyConst<AU>::secondToTime(2 * 1E-12), thermo, kineticModel, forceModel);

            typename MDCell<ScalarType, PosScalarType>::PositionMatrix buffer(forceModel.getNumMolecule(), 3, 0);
            for (size_t i = 0; i < numReplica; ++i)
                buffer += forceModel.makePermanentDipole(rpmd.phaseToCell(i));
            buffer *= reciprocal(ScalarType(numReplica));
            Vector<ScalarType, 3> dipole0{buffer.col(0).asVector().sum(), buffer.col(1).asVector().sum(), buffer.col(2).asVector().sum()};
            for (size_t i = 0; i < corr.getLength(); ++i) {
                buffer = ScalarType(0);
                for (size_t i = 0; i < numReplica; ++i)
                    buffer += forceModel.makePermanentDipole(rpmd.phaseToCell(i));
                buffer *= reciprocal(ScalarType(numReplica));
                Vector<ScalarType, 3> dipole1{buffer.col(0).asVector().sum(), buffer.col(1).asVector().sum(), buffer.col(2).asVector().sum()};
                toNextMean(corr[i], path, dipole0 * dipole1);
                rpmd.nvt_step<ThermostatType, RandomPoolType, KineticModel, ForceModel, ThreadExecutor>(thermo, kineticModel, forceModel);
            }
            rpmd.normalizeCentroid();
            std::cout << path << std::endl;
        }
    }
    ThreadPool::getInstance().shouldExit();

    IRAbsorb<ScalarType> ir(corr, temperatureT, timeStep, double(rpmd.getVolume()), 128, 4);
    auto spectrum = ir.makeSpectrum();
    spectrum *= reciprocal(ScalarType(PhyConst<SI>::bohrRadius * 100));
    auto waveNum = ir.makeWaveNum();
    waveNum *= reciprocal(ScalarType(PhyConst<SI>::bohrRadius * 100));
    {
        std::ofstream fout("spectrum");
        fout << spectrum;
        fout.close();
    }
    {
        std::ofstream fout("waveNum");
        fout << waveNum;
        fout.close();
    }
    {
        std::ofstream fout("corr");
        fout << corr;
        fout.close();
    }
    {
        Poscar<PosScalarType> poscar({rpmd.getLattice(), rpmd.getRingPolymer().makeCentroidPos(), CrystalCell<PosScalarType>::Type::Cartesian}, {1, 8}, {rpmd.getNumParticle() * 2 / 3, rpmd.getNumParticle() / 3});
        std::ofstream fout("H2O.vasp");
        fout << poscar;
        fout.close();
    }
}
/**
 * Reference:
 * [1] S. Habershon, T. E. Markland, and D. E. Manolopoulosa, J. Chem. Phys. 131, 024501(2009)
 */
int main() {
    testMD();
    return 0;
}
