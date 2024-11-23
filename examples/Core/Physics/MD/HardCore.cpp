/*
 * Copyright 2023 Weibo He.
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
#include <QApplication>
#include <QtCharts/QValueAxis>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Random/RandomSeed.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, 1, 1>;
using MDCellType = MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using KineticModel = HardCore<ScalarType, true, 1, RPMDIntegrator::Exact>;
using RandomType = Random<MT19937>;
constexpr double timeStep = 0.001;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 20;
constexpr size_t numMolecular = 20;
constexpr double energy = 20;
constexpr double temperatureT = 2 * energy / numMolecular;

MDCellType makeSystem() {
    MDCellType::LatticeMatrix lattice{latticeSize};

    auto posVec = VectorND<ScalarType>::random_uniform<RandomType>(numMolecular);
    posVec *= ScalarType(latticeSize);
    std::sort(posVec.begin(), posVec.end());
    MDCellType::PositionMatrix pos(numMolecular, 1);
    pos.col(0) = posVec;

    auto massVec = MDCellType::MassVector::random_uniform<RandomType>(numMolecular);
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

int main(int argc, char** argv) {
    MDType rpmd = MDType(makeSystem(), 1, 1, 1, timeStep);
    rpmd.initMomentum<KineticModel, RandomType>();
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 1, 100);
    kineticModel.updateMass(rpmd.getRingPolymer());

    MatrixType record(10000, rpmd.getNumParticle());
    for (size_t i = 0; i < record.getRow(); ++i) {
        ForceModel forceModel{};
        rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, forceModel);
        record.row(i) = rpmd.getRingPolymer().makeBeadPos(0).col(0);
    }
    const VectorType t = VectorType::linspace(0, record.getRow() * timeStep, record.getRow());

    QApplication app(argc, argv);
    constexpr double minX = 0;
    constexpr double maxX = 10.1;
    Plot* plot = new Plot(minX, maxX, 0 - 5, 20 + 5, 2, 5);
    auto& chart = *plot->getChart();
    chart.legend()->setVisible(false);
    plot->getAxisX()->setLabelFormat("%d");
    plot->getAxisY()->setLabelFormat("%d");
    plot->getAxisX()->setTitleText("Time");
    plot->getAxisY()->setTitleText("x");

    for (size_t i = 0; i < rpmd.getNumParticle(); ++i)
        plot->line(t, record.col(i));
    auto& line1 = plot->line(VectorType{minX, maxX}, VectorType{0, 0});
    line1.setColor(Qt::black);
    auto& line2 = plot->line(VectorType{minX, maxX}, VectorType{latticeSize, latticeSize});
    line2.setColor(Qt::black);

    plot->show();
    return QApplication::exec();
}
