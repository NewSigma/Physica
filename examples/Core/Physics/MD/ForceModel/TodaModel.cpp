/*
 * Copyright 2023-2024 Weibo He.
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
#include <QtCharts/QValueAxis>
#include <gperftools/profiler.h>
#include "Physica/Core/Math/Calculus/Integrate/Integrate.h"
#include "Physica/Core/Math/Random/RandomSeed.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/Langevin.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Physics/MD/ForceModel/EmptyForceModel.h"
#include "Physica/Core/Physics/MD/ForceModel/TodaModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using ThermoType = Langevin<ScalarType, 1, 1>;
using KineticModel = HardCore<ScalarType, false, 1, RPMDIntegrator::Exact>;
using ForceModel = TodaModel<ScalarType, true>;
using MDType = RPMD<ScalarType, 1, 1>;
using MDCellType = MDType::MDCellType;
using RandomType = Random<std::mt19937>;
constexpr size_t numReplica = 1;
constexpr double temperatureT = 0.02;
constexpr double thermostatTime = 0.01;
constexpr double timeStep = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize0 = 20;
constexpr unsigned int numMolecular = 20;

MDCellType makeSystem(ScalarType latticeSize) {
    MDCellType::LatticeMatrix lattice{latticeSize};

    MDCellType::PositionMatrix pos(numMolecular, 1);
    for (size_t i = 0; i < numMolecular; ++i)
        pos(i, 0) = latticeSize * ScalarType(i + 1) / ScalarType(numMolecular + 1);

    MDCellType::MassVector massVec(numMolecular);
    for (size_t i = 0; i < numMolecular; ++i) {
        massVec[i] = (i % 2 == 0) ? 1.0 : 2.0;
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

std::pair<ScalarType, ScalarType> calcPress(size_t numSystem, size_t numStep, ScalarType latticeSize) {
    auto rpmd = MDType(makeSystem(latticeSize), numReplica, numReplica, temperatureT, timeStep);
    const ThermoType thermo(temperatureT, thermostatTime, true);
    ForceModel forceModel(1.0);
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 1, 100);
    kineticModel.updateMass(rpmd.getRingPolymer());
    auto& pool = RandomType::getInstance();
    auto& gen = pool.getGen();
    rpmd.initMomentum<KineticModel, decltype(gen)>(gen);

    ScalarType mean = 0, variance = 0;
    for (size_t sys = 0; sys < numSystem; ++sys) {
        ScalarType temp = 0;
        for (size_t i = 0; i < numStep; ++i) {
            rpmd.nvt_step<ThermoType, RandomType, KineticModel, ForceModel, SequentialExecutor>(
                    thermo, pool, kineticModel, forceModel);
            toNextMean(temp, i, rpmd.makeStressClassical<ForceModel, SequentialExecutor>(forceModel)(0, 0));
        }
        toNextVariance(variance, mean, sys, temp);
    }
    return std::make_pair(mean, sqrt(variance));
}

Plot& plotPress(const VectorType& lattices, const VectorType& density) {
    VectorType meanPress(lattices.getLength()), deviaPress(lattices.getLength());
    ThreadExecutor::parallel_for([&meanPress, &deviaPress, &lattices](unsigned int i) {
        auto pair = calcPress(8, 100000, lattices[i]);
        meanPress[i] = pair.first;
        deviaPress[i] = pair.second;
    }, lattices.getLength(), 4).wait();

    QFont font;
    Plot* plot = new Plot(0.6, 1.4, -0.5, 0.5, 0.4, 0.5);
    plot->getChart()->legend()->setVisible(false);
    plot->getAxisX()->setLabelFormat("%.1f");
    plot->getAxisY()->setLabelFormat("%.1f");
    plot->getAxisX()->setTitleText("&rho;");
    plot->getAxisY()->setTitleText("p");
    {
        auto& area = plot->area_center(density, meanPress, deviaPress);
        auto& spline = plot->line(density, meanPress);
        spline.setColor(area.color());

        auto color = area.color();
        color.setAlpha(75);
        area.setColor(color);
    }
    return *plot;
}

Plot& plotDeltaFreeEnergy(const VectorType& lattices, const VectorType& density) {
    VectorType deltaFreeEnergy(lattices.getLength());
    {
        deltaFreeEnergy[lattices.getLength() - 1] = 0;
        ThreadExecutor::parallel_for([&deltaFreeEnergy, &density](unsigned int i) {
            Integrate<IntegrateMethod::Simpson, ScalarType, 1> simpson({{density[i + 1]}, {density[i]}}, 0.0001);
            deltaFreeEnergy[i] = simpson.solve([&](ScalarType rho) -> ScalarType {
                const ScalarType numParticle = numMolecular;
                const ScalarType latticeSize = numParticle / rho;
                return calcPress(1, 10000, latticeSize).first * numParticle / square(rho);
            });
        }, lattices.getLength() - 1, 4).wait();

        for (size_t i = lattices.getLength() - 2; i < lattices.getLength(); --i)
            deltaFreeEnergy[i] += deltaFreeEnergy[i + 1];
        const ScalarType perParticleFactor = reciprocal(ScalarType(numMolecular));
        deltaFreeEnergy *= perParticleFactor;
    }
    QFont font;
    Plot* plot = new Plot(0.6, 0.4, -0.11, 0.01, 0.4, 0.02);
    plot->getChart()->legend()->setVisible(false);
    plot->getAxisX()->setLabelFormat("%.1f");
    plot->getAxisY()->setLabelFormat("%.2f");
    plot->getAxisX()->setTitleText("&rho;");
    plot->getAxisY()->setTitleText("&Delta;f");
    plot->line(density, deltaFreeEnergy);
    return *plot;
}
/**
 * Reference:
 * [1] L. Y. Zhu and J. Wang, Phys. Rev. E 98, 022117 (2018)
 */
int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    QApplication app(argc, argv);
    const auto lattices = VectorType::linspace(0.75 * latticeSize0, 1.5 * latticeSize0, 100);
    const VectorType density = reciprocal(lattices * reciprocal(ScalarType(latticeSize0)));
    auto& plot1 = plotPress(lattices, density);
    auto& plot2 = plotDeltaFreeEnergy(lattices, density);
    plot1.show();
    plot2.show();
    return QApplication::exec();
}
