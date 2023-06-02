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
#include <QApplication>
#include <QtCharts/QValueAxis>
#include "Physica/Gui/Plot/Plot.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "Physica/Utils/Random.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double, false>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, ScalarType, 1, 1>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, ScalarType, 1>;
using KineticModel = HardCore<ScalarType, 1>;
constexpr double timeStep = 0.001;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 20;
constexpr size_t numMolecular = 20;
constexpr double energy = 20;
constexpr double temperatureT = 2 * energy / numMolecular;

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
    for (auto& elem : massVec)
        elem = dist(gen);
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

void scaleVelocity(MDType& rpmd) {
    const ScalarType energy1 = rpmd.getRingPolymer().calcClassicalKinetic();
    auto phase = rpmd.getPhaseMatrix().col(0);
    auto momentum = phase.head(numMolecular);
    momentum *= sqrt(ScalarType(energy) / energy1);
}

int main(int argc, char** argv) {
    std::mt19937::result_type seed;
    Physica::Utils::Random::rdrand(seed);
    std::mt19937 gen(seed);

    MDType rpmd = MDType(makeSystem(gen), 1, 1, 1, timeStep);
    rpmd.initMomentum(gen);
    scaleVelocity(rpmd);
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 1, 100);
    kineticModel.updateMass(rpmd.getRingPolymer());

    MatrixType record(10000, rpmd.getNumParticle());
    for (size_t i = 0; i < record.getRow(); ++i) {
        rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, ForceModel());
        record.row(i) = rpmd.getRingPolymer().makeBeadPos(0);
    }
    const VectorType t = VectorType::linspace(0, record.getRow() * timeStep, record.getRow());

    QApplication app(argc, argv);
    QFont font;
    Plot* plot = new Plot();
    auto& chart = *plot->chart();
    chart.legend()->setVisible(false);
    {
        constexpr double minX = 0;
        constexpr double maxX = 10.1;
        constexpr double minY = 0;
        constexpr double maxY = 20;
        constexpr double deltaX = 2;
        constexpr double deltaY = 5;
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
        axisX->setTitleText("Time");
        axisX->setLabelFormat("%d");
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
        axisY->setTitleText("x");
        axisY->setLabelFormat("%d");
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

        for (size_t i = 0; i < rpmd.getNumParticle(); ++i) {
            auto& spline = plot->line(t, record.col(i));
            spline.attachAxis(axisX);
            spline.attachAxis(axisY);
        }
    }
    plot->show();
    return QApplication::exec();
}
