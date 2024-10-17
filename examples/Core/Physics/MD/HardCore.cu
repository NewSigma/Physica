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
#include <algorithm>
#include <fstream>
#include <iostream>
#include <gperftools/profiler.h>
#include <QApplication>
#include <QtCharts/QValueAxis>
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.cuh"
#include "Physica/Core/Math/Random/RandomPool.h"
#include "Physica/Core/Math/Statistics/ProbDistribution.h"
#include "Physica/Core/Parallel/Executor/ThreadExecutor.h"
#include "Physica/Gui/Plot/Plot.h"
#include "Physica/Core/Utils/BenchmarkHelper.h"
#include <Physica/Core/Utils/Cycler.h>

using namespace Physica::Core;
using namespace Physica::Gui;
constexpr double timeStep = 0.1;
constexpr double collideFactor = 0.005;
const size_t numMolecular = 512;
constexpr double latticeSize = 512;
constexpr double temperatureT = 2;
[[maybe_unused]] constexpr double energy = numMolecular * temperatureT / 2;
constexpr size_t numSample = 5000;
constexpr bool IsComputeMode = true;

using ScalarType = float32;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, 1, 1>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using KineticModel = HardCore<ScalarType, false, 1, RPMDIntegrator::Exact, CUDAExecutor>;
using RandomGenerator = std::mt19937;

MDCellType makeSystem(RandomGenerator& gen) {
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
            auto& gen = RandomPool<RandomGenerator>::getInstance().getGen();
            MDType rpmd = MDType(makeSystem(gen), 1, 1, temperatureT, timeStep);
            KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 100);
            kineticModel.updateMass(rpmd.getRingPolymer());
            rpmd.initMomentum<KineticModel, decltype(gen)>(gen);

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
    Plot* plot = new Plot(0, 3.001, -3, 1, 1, 1);
    auto& chart = *plot->getChart();
    chart.legend()->setVisible(true);
    chart.legend()->setAlignment(Qt::AlignTop);
    chart.legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    {
        auto* axisX = plot->getAxisX();
        auto* axisY = plot->getAxisY();
        axisX->setTitleText("t");
        axisX->setLabelFormat("10<sup>%d</sup>");
        axisY->setTitleText("C<sub>JJ</sub>(t)");
        axisY->setLabelFormat("10<sup>%d</sup>");

        VectorType t = VectorType::linspace(1, ScalarType(record.getRow()), record.getRow() - 1);
        t = ln(t) * reciprocal(ln(ScalarType(10)));
        plot->line(t, lgCorr.tail(1));
    }
    plot->show();
    return QApplication::exec();
}
