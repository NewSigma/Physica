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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Random/RandomSeed.h"
#include "Physica/Core/Math/Statistics/ProbabilityDistributionFunction2D.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Physics/MD/Thermostat/Langevin.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "Physica/Gui/Plot/Plot3D.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, ScalarType, 1, 1>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, ScalarType, 1>;
using ThermostatType = Langevin<ScalarType, ScalarType, 1, 1>;
using KineticModel = HardCore<ScalarType, true, 1, RPMDIntegrator::Exact>;
using RandomPoolType = RandomPool<std::mt19937>;
constexpr double timeStepLambda = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 20;
constexpr size_t numMolecular = 20;
constexpr double temperatureT = 2;
constexpr double thermostatTime = 0.01;
constexpr double unitMassM = 1;
constexpr size_t maxHandleNum = 100;
constexpr size_t numStep = 2000000;

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
        massVec[i] = (i % 2U == 0) ? unitMassM : (unitMassM * 10);
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

int main(int argc, char** argv) {
    const double timeStep = timeStepLambda * (latticeSize / numMolecular) * std::sqrt(unitMassM / temperatureT);
    auto& gen = RandomPoolType::getGen();

    MDType rpmd = MDType(makeSystem(gen), 1, 1, temperatureT, timeStep);
    rpmd.initMomentum(gen);
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 1, maxHandleNum);
    kineticModel.updateMass(rpmd.getRingPolymer());
    ThermostatType thermo(temperatureT, thermostatTime);

    ProbabilityDistributionFunction2D<ScalarType> pdf(-10, 10, -10, 10, 100, 100);
    for (size_t i = 0; i < numStep; ++i) {
        ForceModel forceModel{};
        rpmd.nvt_step<ThermostatType, RandomPoolType, KineticModel, ForceModel, SequentialExecutor>(thermo, kineticModel, forceModel);
        pdf.sample(rpmd.getRingPolymer().asMatrix()(0, 0), rpmd.getRingPolymer().asMatrix()(1, 0));
    }
    const auto grid = pdf.makePosition();
    const auto z = pdf.makeDistribution();

    QApplication app(argc, argv);
    Plot3D* plot3d = new Plot3D();
    auto& surf = plot3d->surf(grid.first, grid.second, z);
    surf.activeTheme()->setType(Q3DTheme::ThemePrimaryColors);
    {
        QLinearGradient gr;
        gr.setColorAt(0.0, Qt::blue);
        gr.setColorAt(0.35, Qt::cyan);
        gr.setColorAt(0.5, Qt::green);
        gr.setColorAt(0.65, Qt::yellow);
        gr.setColorAt(1.0, Qt::red);

        surf.seriesList().at(0)->setBaseGradient(gr);
        surf.seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
    }
    plot3d->show();

    return QApplication::exec();
}
