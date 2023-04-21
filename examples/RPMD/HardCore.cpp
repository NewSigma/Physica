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
using namespace Physica::Core::Parallel;
using namespace Physica::Gui;
using ScalarType = Scalar<Double, false>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, ScalarType, 1>;
using MDCellType = typename MDType::MDCellType;
constexpr double timeStep = 0.001;
constexpr double latticeSize = 100;
constexpr size_t numMolecular = 10;

template<class ScalarType>
class EmptyForceModel final {
public:
    /* Operations */
    template<class Executor>
    [[nodiscard]] Vector<ScalarType> force(const MDCellType& cell) const { return Vector<ScalarType>(cell.getDOF(), 0); }
    template<class Executor>
    [[nodiscard]] Vector<ScalarType> force_short(const MDCellType& cell) const { return force<Executor>(cell); }
    template<class Executor>
    [[nodiscard]] Vector<ScalarType> force_long(const MDCellType& cell) const { return Vector<ScalarType>(cell.getDOF(), 0); }
};

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

int main(int argc, char** argv) {
    using ForceModel = EmptyForceModel<ScalarType>;
    using KineticModel = HardCore<ScalarType>;
    std::mt19937::result_type seed;
    Physica::Utils::Random::rdrand(seed);
    std::mt19937 gen(seed);

    MDType rpmd = MDType(makeSystem(gen), 1, 1, 1, timeStep);
    rpmd.initMomentum(gen);
    KineticModel kineticModel(latticeSize, 0.01);

    MatrixType record(100000, rpmd.getNumParticle());
    for (size_t i = 0; i < record.getRow(); ++i) {
        rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, ForceModel());
        record.row(i) = rpmd.getRingPolymer().makeBeadPos(0);
    }

    QApplication app(argc, argv);
    QFont font;
    Plot* plot = new Plot();
    auto& chart = *plot->chart();
    chart.legend()->setVisible(false);
    {
        constexpr double minX = 0;
        constexpr double maxX = 100000;
        constexpr double minY = 0;
        constexpr double maxY = 100;
        constexpr double deltaX = 20000;
        constexpr double deltaY = 20;
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
        axisX->setTitleText("Step");
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
        axisY->setTitleText("X");
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
            auto& spline = plot->line(record.col(i));
            spline.attachAxis(axisX);
            spline.attachAxis(axisY);
        }
    }
    plot->show();
    return QApplication::exec();
}
