/*
 * Copyright 2023-2025 Weibo He.
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
#include "Physica/Core/Math/Statistics/ProbDistribution2D.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Physics/MD/Thermostat/Langevin.h"
#include "Physica/Core/Parallel/Executor/SeqExecutor.h"
#include "Physica/Gui/Plot/Plot3D.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, 1, 1>;
using MDCellType = MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using ThermoType = Langevin<ScalarType, 1, 1>;
using KineticModel = HardCore<ScalarType, true, 1, RPMDIntegrator::Exact>;
using RandomType = Random<MT19937>;
constexpr double timeStepLambda = 0.01;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 20;
constexpr size_t numMolecular = 20;
constexpr double temperatureT = 2;
constexpr double thermostatTime = 0.01;
constexpr double unitMassM = 1;
constexpr size_t maxHandleNum = 100;
constexpr size_t numStep = 2000000;

MDCellType makeSystem() {
    MDCellType::LatticeMatrix lattice{latticeSize};

    auto posVec = VectorND<ScalarType>::random_uniform<RandomType>(numMolecular);
    posVec *= ScalarType(latticeSize);
    std::sort(posVec.begin(), posVec.end());
    MDCellType::PositionMatrix pos(numMolecular, 1);
    pos.col(0) = posVec;

    MDCellType::MassVector massVec(numMolecular);
    for (size_t i = 0; i < numMolecular; ++i) {
        massVec[i] = (i % 2U == 0) ? unitMassM : (unitMassM * 10);
    }
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

int main(int argc, char** argv) {
    const double timeStep = timeStepLambda * (latticeSize / numMolecular) * std::sqrt(unitMassM / temperatureT);

    MDType rpmd = MDType(makeSystem(), 1, 1, temperatureT, timeStep);
    rpmd.initMomentum<KineticModel, RandomType>();
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 1, maxHandleNum);
    kineticModel.updateMass(rpmd.getRingPolymer());
    ThermoType thermo(temperatureT, thermostatTime, true);

    ProbDistribution2D<ScalarType> pdf(-10, 10, -10, 10, 100, 100);
    for (size_t i = 0; i < numStep; ++i) {
        ForceModel forceModel{};
        rpmd.nvt_step<ThermoType, RandomType, KineticModel, ForceModel, SeqExecutor>(thermo, kineticModel, forceModel);
        pdf.sample(rpmd.getRingPolymer().asMatrix()(0, 0), rpmd.getRingPolymer().asMatrix()(1, 0));
    }
    const auto grid = pdf.makePosition();
    const auto z = pdf.makeDistribution();

    QApplication app(argc, argv);
    Plot3D* plot3d = new Plot3D();
    auto& surf = plot3d->surf(grid.first, grid.second, z);
    surf.setBaseGradient(Plot3D::makeDefaultGrad());
    surf.setColorStyle(Q3DTheme::ColorStyleRangeGradient);
    plot3d->show();
    return QApplication::exec();
}
