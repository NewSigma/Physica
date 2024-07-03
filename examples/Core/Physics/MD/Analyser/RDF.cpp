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
#include <QApplication>
#include "Physica/Core/Physics/MD/ForceModel/Ewald/RandomBatchEwald.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/SolidState/CrystalCell.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/DoubleThermo.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Physics/MD/Analyser/RDF.h"
#include "Physica/Core/IO/VASP/Poscar.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using RandomPoolType = RandomPool<std::mt19937>;
using EwaldType = RandomBatchEwald<ScalarType, RandomPoolType>;
using ForceModel = Q_TIP4P<ScalarType, EwaldType>;
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

template<size_t NumReplica>
RDF<ScalarType> calcRDF(size_t numReplica) {
    using KineticModel = FreeModel<ScalarType, 3, NumReplica, RPMDIntegrator::Exact>;
    using ThermoType = DoubleThermo<KineticModel>;
    auto& pool = RandomPoolType::getInstance();
    auto& gen = pool.getGen();
    auto cell = makeSystem(3, gen);
    ForceModel::sortPosition(cell);
    const ThermoType thermo(temperatureT, thermostatTime);
    RPMD<ScalarType, 3, NumReplica> rpmd(std::move(cell), numReplica, numReplica, temperatureT, timeStep);
    rpmd.template initMomentum<KineticModel, decltype(gen)>(gen);
    ForceModel forceModel(rpmd.phaseToCell(0), pair_cutoff, EwaldType(1000, 100));
    KineticModel kineticModel(temperatureT, numReplica);

    RDF<ScalarType> rdf;
    {
        Physica::Utils::Array<bool> isFromParticle(rpmd.getNumParticle());
        Physica::Utils::Array<bool> isToParticle(rpmd.getNumParticle());
        for (size_t i = 0; i < isFromParticle.getLength(); ++i) {
            const bool isHydrogen = i < isFromParticle.getLength() * 2 / 3;
            isFromParticle[i] = isHydrogen;
            isToParticle[i] = isHydrogen;
        }
        rdf = RDF<ScalarType>(std::move(isFromParticle), std::move(isToParticle), rpmd.getVolume(), PhyConst<AU>::angstormToBohr(0.01), 700);
    }

    ThreadPool::numThreadRequired = 4;
    {
        rpmd.template nvt_step_for<ThermoType, RandomPoolType, KineticModel, ForceModel, ThreadExecutor>(
            PhyConst<AU>::secondToTime(2 * 1E-12), thermo, pool, kineticModel, forceModel);
        for (size_t i = 0; i < 1000; ++i) {
            rpmd.template nvt_step<ThermoType, RandomPoolType, KineticModel, ForceModel, ThreadExecutor>(
                thermo, pool, kineticModel, forceModel);
            for (size_t j = 0; j < numReplica; ++j)
                rdf.sample(rpmd.phaseToCell(j));
        }
    }
    return rdf;
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Plot* plot = new Plot(0, 7, 0, 3, 2, 1);
    auto* legend = plot->getChart()->legend();
    legend->setVisible(true);
    legend->setAlignment(Qt::AlignTop);
    legend->setMarkerShape(QLegend::MarkerShapeFromSeries);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("r/Å");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("g<sub>HH</sub>/Å<sup>-1</sup>");
    axisY->setLabelFormat("%d");
    {
        const auto rdf = calcRDF<Physica::Dynamic>(32);
        const VectorType dists = rdf.makeDists() * ScalarType(PhyConst<AU>::bohrToAngstorm(1));
        const VectorType rdfLine = rdf.makeRDF();
        auto& line = plot->line(dists, rdfLine);
        line.setName("PIMD");
    }
    {
        const auto rdf = calcRDF<1>(1);
        const VectorType dists = rdf.makeDists() * ScalarType(PhyConst<AU>::bohrToAngstorm(1));
        const VectorType rdfLine = rdf.makeRDF();
        auto& line = plot->line(dists, rdfLine);
        line.setName("MD");
    }
    plot->show();
    return QApplication::exec();
}
