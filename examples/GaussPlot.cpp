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
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Parallel/Executor/SequentialExecutor.h"
#include "Physica/Gui/Plot/GaussPlot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, ScalarType, 1, 1>;
using MDCellType = typename MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, ScalarType, 1>;
using KineticModel = HardCore<ScalarType, true, 1, RPMDIntegrator::Exact>;
constexpr double timeStep = 0.001;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 512;
constexpr size_t numMolecular = 512;
constexpr double energy = 512;
constexpr double temperatureT = 2 * energy / numMolecular;
constexpr size_t numStep = 1000000;
constexpr size_t numSystem = 8;

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
    std::mt19937::result_type seed;
    RandomSeed::rdrand(seed);
    std::mt19937 gen(seed);

    MDType rpmd = MDType(makeSystem(gen), 1, 1, temperatureT, timeStep);
    rpmd.initMomentum(gen);
    KineticModel kineticModel(latticeSize, collideFactor, temperatureT, numMolecular, 1, 1000);
    kineticModel.updateMass(rpmd.getRingPolymer());

    const VectorType repMass = reciprocal(rpmd.getMassVec());
    VectorType meanKinetic(numMolecular, 0);
    VectorType varKinetic(numMolecular, 0);
    VectorType temp(numMolecular);
    for (size_t sys = 0; sys < numSystem; ++sys) {
        VectorType buffer(numMolecular, 0);
        for (size_t i = 0; i < numStep; ++i) {
            ForceModel forceModel{};
            rpmd.nve_step<KineticModel, ForceModel, SequentialExecutor>(kineticModel, forceModel);
            auto col = rpmd.getRingPolymer().asMatrix().col(0);
            auto momentum = col.head(numMolecular);
            temp = hadamard(square(momentum), repMass) * ScalarType(0.5);
            toNextMean(buffer, i, temp);
        }

        for (size_t j = 0; j < numMolecular; ++j) {
            toNextVariance(varKinetic[j], meanKinetic[j], sys, buffer[j]);
        }
    }
    const VectorType deviaKinetic = sqrt(varKinetic);

    QApplication app(argc, argv);
    GaussPlot* plot = new GaussPlot(0.25, 0.5, 0.25, 0.25, 2);
    plot->chart()->legend()->hide();
    {
        VectorType moved = meanKinetic - ScalarType(temperatureT * 0.5);
        auto& scatter = plot->scatter(moved, deviaKinetic);
        scatter.setMarkerSize(10);
        scatter.setColor(Qt::red);
    }
    plot->show();
    return QApplication::exec();
}
