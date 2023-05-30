/*
 * Copyright 2023 WeiBo He.
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
#include <gperftools/profiler.h>
#include "Physica/Core/Math/Calculus/Integrate/Integrate.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/Thermostat/Langevin.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Physics/MD/ForceModel/FreeModel.h"
#include "Physica/Core/Physics/MD/ForceModel/TodaModel.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Utils/Random.h"
#include <QApplication>
#include <QtCharts/QValueAxis>
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double, false>;
using PosScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using ThermostatType = Langevin<ScalarType, PosScalarType, 1, 1>;
using KineticModel = HardCore<ScalarType>;
using ForceModel = TodaModel<ScalarType, PosScalarType, 1>;
using MDType = RPMD<ScalarType, PosScalarType, 1, 1>;
using MDCellType = typename MDType::MDCellType;
constexpr size_t numReplica = 1;
constexpr double temperatureT = 0.02;
constexpr double thermostatTime = 0.01;
constexpr double timeStep = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize0 = 20;
constexpr unsigned int numMolecular = 20;

MDCellType makeSystem(ScalarType latticeSize) {
    typename MDCellType::LatticeMatrix lattice{latticeSize};

    typename MDCellType::PositionMatrix pos(numMolecular, 1);
    for (size_t i = 0; i < numMolecular; ++i)
        pos(i, 0) = latticeSize * ScalarType(i + 1) / ScalarType(numMolecular + 1);

    typename MDCellType::MassVector massVec(numMolecular);
    for (size_t i = 0; i < numMolecular; ++i) {
        massVec[i] = (i % 2 == 0) ? 1.0 : 2.0;
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

std::pair<ScalarType, ScalarType> calcPress(size_t numSystem, size_t numStep, ScalarType latticeSize, std::mt19937& gen) {
    auto rpmd = MDType(makeSystem(latticeSize), numReplica, numReplica, temperatureT, timeStep);
    const ThermostatType thermo(temperatureT, thermostatTime);
    const ForceModel forceModel{};
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 100);
    kineticModel.updateMass(rpmd.getRingPolymer());
    rpmd.initMomentum(gen);
    rpmd.updateForce<ForceModel, SequentialExecutor>(forceModel);

    ScalarType mean = 0, variance = 0;
    for (size_t sys = 0; sys < numSystem; ++sys) {
        ScalarType temp = 0;
        for (size_t i = 0; i < numStep; ++i) {
            rpmd.nvt_step<ThermostatType, decltype(gen), KineticModel, ForceModel, SequentialExecutor>(
                    thermo, gen, kineticModel, forceModel);
            toNextMean(temp, i, rpmd.makeStress(forceModel)(0, 0));
        }
        toNextVariance(variance, mean, sys, temp);
    }
    return std::make_pair(mean, sqrt(variance));
}

void plotPress(const VectorType& lattices, const VectorType& density) {
    VectorType meanPress(lattices.getLength()), deviaPress(lattices.getLength());
    ThreadExecutor::parallel_for([&meanPress, &deviaPress, &lattices](unsigned int i) {
        std::mt19937::result_type seed;
        Physica::Utils::Random::rdrand(seed);
        std::mt19937 gen(seed);
        auto pair = calcPress(8, 100000, lattices[i], gen);
        meanPress[i] = pair.first;
        deviaPress[i] = pair.second;
    }, lattices.getLength(), 4).wait();

    QFont font;
    Plot* plot = new Plot();
    auto& chart = *plot->chart();
    chart.legend()->setVisible(false);
    {
        constexpr double minX = 0.6;
        constexpr double maxX = 1.4;
        constexpr double minY = -0.5;
        constexpr double maxY = 0.5;
        constexpr double deltaX = 0.4;
        constexpr double deltaY = 0.5;
        QValueAxis* axisX = new QValueAxis();
        font = axisX->labelsFont();
        font.setPointSize(15);
        axisX->setTickAnchor(0);
        axisX->setTickInterval(deltaX);
        axisX->setTickType(QValueAxis::TicksDynamic);
        axisX->setMinorGridLineVisible(false);
        axisX->setMinorTickCount(4);
        axisX->setLinePenColor(Qt::black);
        axisX->setGridLineVisible(false);
        axisX->setLabelsFont(font);
        axisX->setRange(minX, maxX);
        axisX->setTitleText("&rho;");
        axisX->setLabelFormat("%.1f");
        axisX->setTitleFont(font);
        QValueAxis* axisY = new QValueAxis();
        axisY->setTickAnchor(0);
        axisY->setTickInterval(deltaY);
        axisY->setTickType(QValueAxis::TicksDynamic);
        axisY->setMinorGridLineVisible(false);
        axisY->setMinorTickCount(4);
        axisY->setLinePenColor(Qt::black);
        axisY->setGridLineVisible(false);
        axisY->setMinorGridLineVisible(false);
        axisY->setLabelsFont(font);
        axisY->setRange(minY, maxY);
        axisY->setTitleText("p");
        axisY->setLabelFormat("%.1f");
        axisY->setTitleFont(font);
        QValueAxis* axisTop = new QValueAxis();
        axisTop->setTickAnchor(0);
        axisTop->setTickInterval(deltaX);
        axisTop->setTickType(QValueAxis::TicksDynamic);
        axisTop->setLabelsVisible(false);
        axisTop->setMinorTickCount(4);
        axisTop->setGridLineVisible(false);
        axisTop->setMinorGridLineVisible(false);
        axisTop->setRange(minX, maxX);
        axisTop->setLinePenColor(Qt::black);
        QValueAxis* axisRight = new QValueAxis();
        axisRight->setTickAnchor(0);
        axisRight->setTickInterval(deltaY);
        axisRight->setTickType(QValueAxis::TicksDynamic);
        axisRight->setLabelsVisible(false);
        axisRight->setGridLineVisible(false);
        axisRight->setMinorGridLineVisible(false);
        axisRight->setMinorTickCount(4);
        axisRight->setRange(minY, maxY);
        axisRight->setLinePenColor(Qt::black);

        chart.addAxis(axisX, Qt::AlignBottom);
        chart.addAxis(axisY, Qt::AlignLeft);
        chart.addAxis(axisTop, Qt::AlignTop);
        chart.addAxis(axisRight, Qt::AlignRight);

        {
            auto& area = plot->area_center(density, meanPress, deviaPress);
            area.attachAxis(axisX);
            area.attachAxis(axisY);

            auto& spline = plot->line(density, meanPress);
            spline.setColor(area.color());
            spline.attachAxis(axisX);
            spline.attachAxis(axisY);

            auto color = area.color();
            color.setAlpha(75);
            area.setColor(color);
        }
    }
    plot->show();
}

void plotDeltaFreeEnergy(const VectorType& lattices, const VectorType& density) {
    VectorType deltaFreeEnergy(lattices.getLength());
    {
        deltaFreeEnergy[lattices.getLength() - 1] = 0;
        ThreadExecutor::parallel_for([&deltaFreeEnergy, &density](unsigned int i) {
            std::mt19937::result_type seed;
            Physica::Utils::Random::rdrand(seed);
            std::mt19937 gen(seed);

            Integrate<Simpson, ScalarType, 1> simpson({{density[i + 1]}, {density[i]}}, 0.0001);
            deltaFreeEnergy[i] = simpson.solve([&](ScalarType rho) -> ScalarType {
                const ScalarType numParticle = numMolecular;
                const ScalarType latticeSize = numParticle / rho;
                return calcPress(1, 10000, latticeSize, gen).first * numParticle / square(rho);
            });
        }, lattices.getLength() - 1, 4).wait();

        for (size_t i = lattices.getLength() - 2; i < lattices.getLength(); --i)
            deltaFreeEnergy[i] += deltaFreeEnergy[i + 1];
        const ScalarType perParticleFactor = reciprocal(ScalarType(numMolecular));
        deltaFreeEnergy *= perParticleFactor;
    }
    QFont font;
    Plot* plot = new Plot();
    auto& chart = *plot->chart();
    chart.legend()->setVisible(false);
    {
        constexpr double minX = 0.6;
        constexpr double maxX = 1.4;
        constexpr double minY = -0.11;
        constexpr double maxY = 0.01;
        constexpr double deltaX = 0.4;
        constexpr double deltaY = 0.02;
        QValueAxis* axisX = new QValueAxis();
        font = axisX->labelsFont();
        font.setPointSize(15);
        axisX->setTickAnchor(0);
        axisX->setTickInterval(deltaX);
        axisX->setTickType(QValueAxis::TicksDynamic);
        axisX->setMinorGridLineVisible(false);
        axisX->setMinorTickCount(4);
        axisX->setLinePenColor(Qt::black);
        axisX->setGridLineVisible(false);
        axisX->setLabelsFont(font);
        axisX->setRange(minX, maxX);
        axisX->setTitleText("&rho;");
        axisX->setLabelFormat("%.1f");
        axisX->setTitleFont(font);
        QValueAxis* axisY = new QValueAxis();
        axisY->setTickAnchor(0);
        axisY->setTickInterval(deltaY);
        axisY->setTickType(QValueAxis::TicksDynamic);
        axisY->setMinorGridLineVisible(false);
        axisY->setMinorTickCount(4);
        axisY->setLinePenColor(Qt::black);
        axisY->setGridLineVisible(false);
        axisY->setMinorGridLineVisible(false);
        axisY->setLabelsFont(font);
        axisY->setRange(minY, maxY);
        axisY->setTitleText("&Delta;f");
        axisY->setLabelFormat("%.2f");
        axisY->setTitleFont(font);
        QValueAxis* axisTop = new QValueAxis();
        axisTop->setTickAnchor(0);
        axisTop->setTickInterval(deltaX);
        axisTop->setTickType(QValueAxis::TicksDynamic);
        axisTop->setLabelsVisible(false);
        axisTop->setMinorTickCount(4);
        axisTop->setGridLineVisible(false);
        axisTop->setMinorGridLineVisible(false);
        axisTop->setRange(minX, maxX);
        axisTop->setLinePenColor(Qt::black);
        QValueAxis* axisRight = new QValueAxis();
        axisRight->setTickAnchor(0);
        axisRight->setTickInterval(deltaY);
        axisRight->setTickType(QValueAxis::TicksDynamic);
        axisRight->setLabelsVisible(false);
        axisRight->setGridLineVisible(false);
        axisRight->setMinorGridLineVisible(false);
        axisRight->setMinorTickCount(4);
        axisRight->setRange(minY, maxY);
        axisRight->setLinePenColor(Qt::black);

        chart.addAxis(axisX, Qt::AlignBottom);
        chart.addAxis(axisY, Qt::AlignLeft);
        chart.addAxis(axisTop, Qt::AlignTop);
        chart.addAxis(axisRight, Qt::AlignRight);

        auto& spline = plot->line(density, deltaFreeEnergy);
        spline.attachAxis(axisX);
        spline.attachAxis(axisY);
    }
    plot->show();
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
    plotPress(lattices, density);
    plotDeltaFreeEnergy(lattices, density);
    return QApplication::exec();
}
