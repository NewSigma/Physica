/*
 * Copyright 2024 Weibo He.
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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using namespace Physica::Gui;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using MatrixType = DenseMatrix<ScalarType>;
constexpr int NumAve = 100;
constexpr int NumStep = 41;
constexpr int NumSite = 4;
constexpr int NumState = 256;

VectorType readTPQ() {
    VectorType result(NumStep);
    MatrixType buffer(NumStep, 6);
    for (int i = 0; i < NumAve; ++i) {
        char path[32];
        sprintf(path, "SS_rand%d.dat", i);
        std::ifstream fin(path);
        if (!fin)
            exit(EXIT_FAILURE);
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        for (size_t r = 0; r < buffer.getRow(); ++r)
            for (size_t c = 0; c < buffer.getCol(); ++c)
                fin >> buffer(r, c);
        toNextMean(result, i, buffer.col(4));
    }
    result *= reciprocal(ScalarType(NumSite));
    return result;
}

VectorType readFullDiag(const VectorType& betas) {
    VectorType energys(NumState), numParticle(NumState);
    {
        MatrixType buffer(NumState, 5);
        std::ifstream fin("zvo_phys.dat");
        if (!fin)
            exit(EXIT_FAILURE);
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        for (size_t r = 0; r < buffer.getRow(); ++r)
            for (size_t c = 0; c < buffer.getCol(); ++c)
                fin >> buffer(r, c);
        energys = buffer.col(0);
        numParticle = buffer.col(1);
    }
    VectorType result(betas.getLength());
    for (size_t i = 0; i < betas.getLength(); ++i) {
        const VectorType weights = exp(energys * (-betas[i]));
        result[i] = numParticle * weights / weights.sum();
    }
    result *= reciprocal(ScalarType(NumSite));
    return result;
}

int main(int argc, char** argv) {
    const auto betas = VectorType::linspace(0, 4, 41);
    const auto rho1 = readTPQ();
    const auto rho2 = readFullDiag(betas);

    QApplication app(argc, argv);
    Plot* plot = new Plot(0, 4, 0.5, 1.05, 1, 0.2);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("&beta;");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("&rho;");
    axisY->setLabelFormat("%.1f");
    plot->line(betas, rho1).setName("TPQ");
    plot->line(betas, rho2).setName("FullDiag");
    plot->show();
    return QApplication::exec();
}
