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
#include "Physica/Core/Physics/MD/KineticModel/HardCore.h"
#include "Physica/Core/Physics/MD/RPMD.h"

import Physica.Gui.GaussPlot;

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
using MDType = RPMD<ScalarType, 1, 1>;
using MDCellType = MDType::MDCellType;
using ForceModel = EmptyForceModel<ScalarType, 1>;
using KineticModel = HardCore<ScalarType, true, 1, RPMDIntegrator::Exact>;
using RandomSource = Random<>;
constexpr double timeStep = 0.001;
constexpr double collideFactor = 0.01;
constexpr double latticeSize = 512;
constexpr size_t numMolecular = 512;
constexpr double energy = 512;
constexpr double temperatureT = 2 * energy / numMolecular;
constexpr size_t numStep = 1000000;
constexpr size_t numSystem = 8;

MDCellType makeSystem() {
    MDCellType::LatticeMatrix lattice{latticeSize};
    auto posVec = VectorND<ScalarType>::random_uniform<RandomSource>(numMolecular);
    posVec *= ScalarType(latticeSize);
    std::sort(posVec.begin(), posVec.end());
    MDCellType::PositionMatrix pos(numMolecular, 1);
    pos.col(0) = posVec;

    auto massVec = MDCellType::MassVector::random_uniform<RandomSource>(numMolecular);
    return MDCellType(std::move(lattice), std::move(pos), std::move(massVec));
}

int main(int argc, char** argv) {
    MDType rpmd = MDType(makeSystem(), 1, 1, temperatureT, timeStep);
    rpmd.initMomentum<KineticModel, RandomSource>();
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
            rpmd.nve_step<Sequential>(kineticModel, forceModel);
            auto col = rpmd.getRingPolymer().asMatrix().col(0);
            auto momentum = col.head(numMolecular);
            temp = hadamard(square(momentum), repMass) * ScalarType(0.5);
            buffer.toNextMean(i, temp);
        }

        for (size_t j = 0; j < numMolecular; ++j) {
            varKinetic[j].toNextVariance(meanKinetic[j], sys, buffer[j]);
        }
    }
    const VectorType deviaKinetic = sqrt(varKinetic);

    QApplication app(argc, argv);
    GaussPlot* plot = new GaussPlot(0.25, 0.5, 0.25, 0.25, 2);
    plot->getChart()->legend()->hide();
    {
        VectorType moved = meanKinetic - ScalarType(temperatureT * 0.5);
        auto& scatter = plot->scatter(moved, deviaKinetic);
        scatter.setMarkerSize(10);
        scatter.setColor(Qt::red);
    }
    plot->show();
    return QApplication::exec();
}
