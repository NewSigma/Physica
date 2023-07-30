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
#include <algorithm>
#include <fstream>
#include <iostream>
#include <gperftools/profiler.h>
#include <QApplication>
#include <QtCharts/QValueAxis>
#include "Physica/Gui/Plot/Plot.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.cuh"
#include "Physica/Utils/Random.h"
#include "Physica/Utils/BenchmarkHelper.h"
#include "Physica/Utils/Cycler.h"
#include "Physica/Utils/CUDA/DeviceProp.cuh"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Core/Math/Statistics/ProbabilityDistributionFunction.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using namespace Physica::Utils;
constexpr double timeStep = 0.1;
constexpr double collideFactor = 0.005;
const size_t numMolecular = 512;
constexpr double latticeSize = 512;
constexpr double temperatureT = 2;
[[maybe_unused]] constexpr double energy = numMolecular * temperatureT / 2;
constexpr size_t numSample = 5000;
constexpr bool IsComputeMode = true;

using ScalarType = Scalar<Float>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, ScalarType, 1, 1>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, ScalarType, 1>;
using KineticModel = HardCore<ScalarType, false, 1, CudaExecutor>;

MDCellType makeSystem(std::mt19937& gen) {
    typename MDCellType::LatticeMatrix lattice{latticeSize};

    std::uniform_real_distribution dist{};
    Vector<ScalarType> posVec(numMolecular);
    for (auto& elem : posVec)
        elem = dist(gen) * latticeSize;
    std::sort(posVec.begin(), posVec.end());
    typename MDCellType::PositionMatrix pos(numMolecular, 1);
    pos.col(0) = posVec;

    typename MDCellType::MassVector massVec(numMolecular);
    for (size_t i = 0; i < numMolecular; ++i) {
        massVec[i] = i % 2 == 0 ? 3.0 : 1.0;
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

ScalarType calcThermoFlux(MDType& rpmd) {
    ScalarType flux = 0;
    auto col = rpmd.getPhaseMatrix().col(0);
    const auto& massVec = rpmd.getMassVec();
    for (size_t i = 0; i < numMolecular; ++i)
        flux += col[i] * square(col[i]) / square(massVec[i]);
    return flux;
}

int main(int argc, char** argv) {
    MatrixType record(1000, 8);
    VectorType mean(record.getRow()), devia(record.getRow());
    if constexpr (IsComputeMode) {
        ThreadPool::numThreadRequired = record.getColumn();
        ThreadExecutor::parallel_for([&record](unsigned int sys) {
            std::mt19937::result_type seed;
            Physica::Utils::Random::rdrand(seed);
            std::mt19937 gen(seed);

            MDType rpmd = MDType(makeSystem(gen), 1, 1, temperatureT, timeStep);
            KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 100);
            kineticModel.updateMass(rpmd.getRingPolymer());
            rpmd.initMomentum(gen);

            Vector<ScalarType> mean(record.getRow(), 0);
            for (size_t sample = 0; sample < numSample; ++sample) {
                const ScalarType flux0 = calcThermoFlux(rpmd);
                kineticModel.nve_step_for(1.0, rpmd.getRingPolymer(), timeStep);
                for (size_t j = 0; j < mean.getLength(); ++j) {
                    const ScalarType flux = calcThermoFlux(rpmd);
                    toNextMean(mean[j], sample, flux0 * flux);
                    kineticModel.do_nve_step_for(1.0, timeStep);
                    kineticModel.post_nve_step(rpmd.getRingPolymer());
                    rpmd.getRingPolymer().removeDrift();
                    kineticModel.updateMomentum(rpmd.getRingPolymer());
                }
            }
            record.asArray()[sys] = std::move(mean);
        }, record.getColumn(), ThreadPool::numThreadRequired).wait();

        for (size_t i = 0; i < mean.getLength(); ++i) {
            mean[i] = Physica::Core::mean(record.row(i));
            devia[i] = Physica::Core::deviation(record.row(i));
        }
        const ScalarType factor = reciprocal(mean[0]);
        mean *= factor;
        devia *= factor;
    }

    if constexpr (!IsComputeMode) {
        std::ifstream fin("data");
        fin >> mean >> devia;
    }
    else {
        std::ofstream fout("data");
        fout << mean << devia;
    }
    const VectorType lgCorr = ln(abs(mean)) * reciprocal(ln(ScalarType(10)));

    QApplication app(argc, argv);
    QFont font;
    Plot* plot = new Plot();
    auto& chart = *plot->chart();
    chart.legend()->setVisible(true);
    chart.legend()->setAlignment(Qt::AlignTop);
    chart.legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    {
        double minX = 0;
        double maxX = 3.001;
        constexpr double minY = -3;
        constexpr double maxY = 1;
        constexpr double deltaX = 1;
        constexpr double deltaY = 1;
        QValueAxis* axisX = new QValueAxis();
        font = axisX->labelsFont();
        font.setPointSize(15);
        axisX->setTickAnchor(0);
        axisX->setTickInterval(deltaX);
        axisX->setTickType(QValueAxis::TicksDynamic);
        axisX->setMinorGridLineVisible(false);
        axisX->setLinePenColor(Qt::black);
        axisX->setGridLineVisible(false);
        axisX->setLabelsFont(font);
        axisX->setRange(minX, maxX);
        axisX->setTitleText("t");
        axisX->setLabelFormat("10<sup>%d</sup>");
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
        axisY->setTitleText("C<sub>JJ</sub>(t)");
        axisY->setLabelFormat("10<sup>%d</sup>");
        axisY->setTitleFont(font);
        QValueAxis* axisTop = new QValueAxis();
        axisTop->setTickAnchor(0);
        axisTop->setTickInterval(deltaX);
        axisTop->setTickType(QValueAxis::TicksDynamic);
        axisTop->setLabelsVisible(false);
        axisTop->setGridLineVisible(false);
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

        VectorType t = VectorType::linspace(1, ScalarType(record.getRow()), record.getRow() - 1);
        t = ln(t) * reciprocal(ln(ScalarType(10)));
        {
            auto& spline = plot->line(t, lgCorr.tail(1));
            spline.attachAxis(axisX);
            spline.attachAxis(axisY);
        }
        /*{
            auto& area = plot->area_center(t, mean, devia);
            area.attachAxis(axisX);
            area.attachAxis(axisY);

            auto& spline = plot->line(t, mean.tail(1));
            spline.setColor(area.color());
            spline.attachAxis(axisX);
            spline.attachAxis(axisY);

            auto color = area.color();
            color.setAlpha(75);
            area.setColor(color);
            spline.setName("Simulation");
        }*/
    }
    plot->show();
    return QApplication::exec();
}
