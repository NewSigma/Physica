/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <format>
#include <fstream>
#include <QApplication>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using T = float64;
using MatrixType = DenseMatrix<T>;
constexpr int NumSystem = 8;
constexpr int NumSample = 100;
constexpr int MaxBeta = 16;
constexpr int NumStep = MaxBeta * 10 + 1;
constexpr int NumSite = 4;
constexpr int NumState = 256;

namespace {
    VectorND<T> readTPQ(int sys) {
        VectorND<T> result(NumStep);
        MatrixType buffer(NumStep, 6);
        for (int i = 0; i < NumSample; ++i) {
            std::ifstream fin(std::format("./{}/output/SS_rand{}.dat", sys, i));
            if (!fin)
                exit(EXIT_FAILURE);
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (size_t r = 0; r < buffer.getRow(); ++r)
                for (size_t c = 0; c < buffer.getCol(); ++c)
                    fin >> buffer[r, c];
            result.toNextMean(i, buffer.col(4));
        }
        result *= reciprocal(T(NumSite));
        return result;
    }

    VectorND<T> readFullDiag(const VectorND<T>& betas) {
        VectorND<T> energys(NumState);
        VectorND<T> numParticle(NumState);
        {
            MatrixType buffer(NumState, 5);
            std::ifstream fin("zvo_phys.dat");
            if (!fin)
                exit(EXIT_FAILURE);
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            for (size_t r = 0; r < buffer.getRow(); ++r)
                for (size_t c = 0; c < buffer.getCol(); ++c)
                    fin >> buffer[r, c];
            energys = buffer.col(0);
            numParticle = buffer.col(1);
        }
        VectorND<T> result(betas.getLength());
        for (size_t i = 0; i < betas.getLength(); ++i) {
            const VectorND<T> weights = exp(energys * (-betas[i]));
            result[i] = numParticle * weights / weights.sum();
        }
        result *= reciprocal(T(NumSite));
        return result;
    }
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Plot* plot = new Plot(0, MaxBeta, 0.5, 1.05, 1, 0.2);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("&beta;");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("&rho;");
    axisY->setLabelFormat("%.1f");

    const auto betas = VectorND<T>::linspace(0, MaxBeta, NumStep);
    {
        VectorND<T> means(betas.getLength(), 0);
        VectorND<T> devias(betas.getLength(), 0);
        for (int i = 0; i < NumSystem; ++i)
            devias.toNextVariance(means, i, readTPQ(i));
        devias = sqrt(devias) * reciprocal(sqrt(T(NumSystem)));

        auto& area = plot->area_center(betas, means, devias * T(2));
        auto color = area.color();
        auto color1 = color;
        color1.setAlpha(75);
        area.setColor(color1);
        auto& l = plot->line(betas, means);
        l.setColor(color);
        l.setName("TPQ");
    }
    plot->line(betas, readFullDiag(betas)).setName("FullDiag");
    plot->show();
    return QApplication::exec();
}
