/*
 * Copyright 2024 WeiBo He.
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
#include <QApplication>
#include "Physica/Core/Math/Random/RandomPool.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/JacobiDavidson.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Physica/Core/Physics/ManyBody/Hubbard1D.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using RandomPoolType = RandomPool<std::mt19937, 10000>;
constexpr unsigned int NumSite = 8;
constexpr unsigned int NumSpinUp = 1;
constexpr unsigned int NumSpinDown = 1;
constexpr double HoppingT = 1.0;
constexpr double RepelU = 2;

int main(int argc, char** argv) {
    Hubbard1D<ScalarType> model({{NumSite}, 1}, HoppingT, RepelU, NumSpinUp, NumSpinDown);
    const size_t numState = model.getNumState();

    QApplication app(argc, argv);
    Plot* plot = new Plot(-0.01, 1.01, -4, 5, 0.2, 2);
    auto& legend = *plot->getChart()->legend();
    legend.setAlignment(Qt::AlignTop);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("N/D");
    axisX->setLabelFormat("%.1f");
    axisY->setTitleText("E");
    axisY->setLabelFormat("%d");
    const auto x = VectorType::linspace(0, 1, numState);
    {
        JacobiDavidson<ScalarType> jd(numState, 48);
        jd.compute(model, VectorType::random_uniform(numState, RandomPoolType::getInstance().getGen()));
        jd.sort();

        auto& scatter = plot->scatter(x.head(jd.getNumRequired()), jd.getEigenvalues());
        scatter.setMarkerSize(10);
        scatter.setName("Sparse");
    }
    {
        using MatrixType = DenseMatrix<ScalarType>;
        const MatrixType h = model;
        SymmEigenSolver<MatrixType> eig(numState);
        eig.compute(h, false);
        eig.sort();

        auto& scatter = plot->scatter(x, eig.getEigenvalues());
        scatter.setMarkerSize(5);
        scatter.setName("Dense");
        scatter.setMarkerShape(QScatterSeries::MarkerShape::MarkerShapeCross);
        auto pen = scatter.pen();
        pen.setColor(Qt::red);
        scatter.setPen(pen);
        scatter.setColor(Qt::red);
    }
    plot->show();
    return QApplication::exec();
}
