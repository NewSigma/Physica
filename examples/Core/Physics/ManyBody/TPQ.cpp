/*
 * Copyright 2024-2026 Weibo He.
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
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/FermiRepr.h"
#include "Physica/Core/Physics/ManyBody/TPQ.h"
#include "Physica/Core/Utils/Container/SymmArray.h"

import Physica.Gui.Plot;

using namespace Physica;
using T = float64;
using RandomSource = Random<>;
constexpr unsigned int NumSite = 4;
constexpr double HoppingT = 1;
constexpr double RepelU = 8;
constexpr double MaxBeta = 16;
constexpr int NumSystem = 8;
constexpr int NumSample = 8192;
constexpr int NumBeta = MaxBeta * 10 + 1;

constexpr int PlotRange = 4;
constexpr int PlotPoint = PlotRange * 10 + 1;

namespace {
    template<class ReprType>
    VectorND<T> calcPartition(ReprType repr_, const VectorND<T>& betas) {
        using Hamilton = HubbardMatrix<T, ReprType>;
        SquareLattice<1> lattice({NumSite}, 1);
        const Hamilton hamilton(HoppingT, RepelU, lattice, std::move(repr_));
        VectorND<T> result(NumBeta);
        const T deltaBeta = betas[1] - betas[0];
        auto psi = TPQ<T>(hamilton.getNumState());
        psi.pre_nvt_step(hamilton, deltaBeta);
        for (unsigned int i = 0; i < NumSample; ++i) {
            psi.template random_normal<RandomSource>();
            for (unsigned int j = 0; j < NumBeta; ++j) {
                if (j != 0)
                    psi.nvt_step(hamilton, deltaBeta);
                result[j].toNextMean(i, exp(psi.lnPartitionZ()));
            }
        }
        return result;
    }

    SymmArray<VectorND<T>> makePartitionMatrix(const VectorND<T>& betas) {
        SymmArray<VectorND<T>> result(NumSite + 1);
        for (int numSpinUp = 0; numSpinUp <= NumSite; ++numSpinUp) {
            for (int numSpinDown = numSpinUp; numSpinDown <= NumSite; ++numSpinDown) {
                const bool useInversionSymm = numSpinUp == numSpinDown;
                VectorND<T> elem;
                if (useInversionSymm) {
                    using ReprType = FermiRepr<1, NumSite, true>;
                    elem = calcPartition(ReprType(numSpinUp, numSpinDown), betas);
                }
                else {
                    using ReprType = FermiRepr<1, NumSite, false>;
                    elem = calcPartition(ReprType(numSpinUp, numSpinDown), betas);
                }
                result[numSpinUp, numSpinDown] = std::move(elem);
            }
        }
        return result;
    }

    VectorND<T> makeDensity(const VectorND<T>& betas) {
        const auto partitionMatrix = makePartitionMatrix(betas);
        VectorND<T> result(NumBeta);
        for (size_t i = 0; i < NumBeta; ++i) {
            T sumPartition = 0;
            T numParticle = 0;
            for (unsigned int numSpinUp = 0; numSpinUp <= NumSite; ++numSpinUp) {
                for (unsigned int numSpinDown = numSpinUp; numSpinDown <= NumSite; ++numSpinDown) {
                    const int factor = numSpinUp == numSpinDown ? 1 : 2;
                    const T n = (numSpinUp + numSpinDown) * factor;
                    sumPartition += partitionMatrix[numSpinUp, numSpinDown][i] * factor;
                    numParticle += partitionMatrix[numSpinUp, numSpinDown][i] * n;
                }
            }
            result[i] = numParticle / (sumPartition * T(NumSite));
        }
        return result;
    }

    std::pair<VectorND<T>, VectorND<T>> makeDensity2(const VectorND<T>& betas) {
        VectorND<T> means(betas.getLength(), 0);
        VectorND<T> devias(betas.getLength(), 0);
        for (int i = 0; i < NumSystem; ++i)
            devias.toNextVariance(means, i, makeDensity(betas));
        devias = sqrt(devias) * reciprocal(sqrt(T(NumSystem)));
        return std::make_pair(std::move(means), std::move(devias));
    }

    void plotTPQ() {
        Plot* plot = new Plot(0, MaxBeta, 0.5, 1.05, 1, 0.2);
        auto* axisX = plot->getAxisX();
        auto* axisY = plot->getAxisY();
        axisX->setTitleText("&beta;");
        axisX->setLabelFormat("%d");
        axisY->setTitleText("&rho;");
        axisY->setLabelFormat("%.1f");

        const auto betas = VectorND<T>::linspace(0, MaxBeta, NumBeta);
        auto [means, devias] = makeDensity2(betas);
        auto& area = plot->area_center(betas, means, devias * T(2));
        auto color = area.color();
        auto color1 = color;
        color1.setAlpha(75);
        area.setColor(color1);
        plot->line(betas, means).setColor(color);
        plot->show();
        QApplication::exec();
    }

    void plotAll() {
        Plot* plot = new Plot(0, PlotRange, 0.5, 1.05, 1, 0.2);
        plot->getLegend().setAlignment(Qt::AlignRight);
        plot->getLegend().show();
        auto* axisX = plot->getAxisX();
        auto* axisY = plot->getAxisY();
        axisX->setTitleText("&beta;");
        axisX->setLabelFormat("%d");
        axisY->setTitleText("&rho;");
        axisY->setLabelFormat("%.1f");

        const auto betas = VectorND<T>::linspace(0, PlotRange, PlotPoint);
        auto h5f = H5File::open("TPQ.h5", H5File::ReadOnly);
        QPen pen1, pen2;
        {
            VectorND<T> tpq_hphi;
            tpq_hphi.read(h5f, "TPQ_HPhi");

            auto& scatter_tpq = plot->scatter(betas, tpq_hphi.head(PlotPoint));
            scatter_tpq.setName("TPQ(HPhi)");
            scatter_tpq.setMarkerSize(10);
            pen1 = scatter_tpq.pen();
            pen1.setColor(scatter_tpq.color());
            scatter_tpq.setPen(pen1);
            scatter_tpq.setColor(Qt::transparent);
        }
        {
            VectorND<T> fulldiag_hphi;
            fulldiag_hphi.read(h5f, "FullDiag_HPhi");

            auto& scatter_fd = plot->scatter(betas, fulldiag_hphi.head(PlotPoint));
            scatter_fd.setName("Full-ED(HPhi)");
            scatter_fd.setMarkerSize(10);
            pen2 = scatter_fd.pen();
            pen2.setColor(scatter_fd.color());
            scatter_fd.setPen(pen2);
            scatter_fd.setColor(Qt::transparent);
        }
        {
            VectorND<T> tpq;
            tpq.read(h5f, "TPQ");

            auto& line_tpq = plot->line(betas, tpq.head(PlotPoint));
            line_tpq.setName("TPQ");
            line_tpq.setColor(pen1.color());
        }
        {
            VectorND<T> fulldiag;
            fulldiag.read(h5f, "FullDiag");

            auto& line_fd = plot->line(betas, fulldiag.head(PlotPoint));
            line_fd.setName("Full-ED");
            auto pen = line_fd.pen();
            pen.setColor(pen2.color());
            pen.setStyle(Qt::DashLine);
            QList<qreal> dashes;
            dashes << 20 << 20;
            pen.setDashPattern(dashes);
            line_fd.setPen(pen);
        }
        plot->show();
        QApplication::exec();
    }
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    plotTPQ();
    plotAll();
    return 0;
}
