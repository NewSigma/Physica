/*
 * Copyright 2024-2025 Weibo He.
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
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/FermiRepr.h"
#include "Physica/Core/Physics/ManyBody/TPQ.h"
#include "Physica/Core/Utils/Container/SymmArray.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using RandomSource = Random<MT19937, 10000>;
constexpr unsigned int NumSite = 4;
constexpr double HoppingT = 1;
constexpr double RepelU = 8;
constexpr unsigned int NumSample = 1000;
constexpr unsigned int NumBeta = 41;

template<class ReprType>
VectorType calcPartition(ReprType repr_, const VectorType& betas) {
    using Hamilton = HubbardMatrix<ScalarType, ReprType>;
    LatticeModel<1> lattice({NumSite}, 1);
    Hubbard<ScalarType, 1> hubbard(lattice, HoppingT, RepelU);
    const Hamilton hamilton(hubbard, std::move(repr_));
    VectorType result(NumBeta);
    const ScalarType deltaBeta = betas[1] - betas[0];
    auto psi = TPQ<ScalarType>(hamilton.getNumState());
    psi.pre_nvt_step(hamilton, deltaBeta);
    for (unsigned int i = 0; i < NumSample; ++i) {
        psi.template random_normal<RandomSource>(1);
        for (unsigned int j = 0; j < NumBeta; ++j) {
            if (j != 0)
                psi.nvt_step(hamilton, deltaBeta);
            toNextMean(result[j], i, psi.calcPartitionXi());
        }
    }
    return result;
}

SymmArray<VectorType> makePartitionMatrix(const VectorType& betas) {
    SymmArray<VectorType> result(NumSite + 1);
    for (unsigned int numSpinUp = 0; numSpinUp <= NumSite; ++numSpinUp) {
        for (unsigned int numSpinDown = numSpinUp; numSpinDown <= NumSite; ++numSpinDown) {
            const bool useInversionSymm = numSpinUp == numSpinDown;
            VectorType elem;
            if (useInversionSymm) {
                using ReprType = FermiRepr<1, NumSite, true>;
                elem = calcPartition(ReprType(numSpinUp, numSpinDown), betas);
            }
            else {
                using ReprType = FermiRepr<1, NumSite, false>;
                elem = calcPartition(ReprType(numSpinUp, numSpinDown), betas);
            }
            result(numSpinUp, numSpinDown) = std::move(elem);
        }
    }
    return result;
}

VectorType makeDensity(const VectorType& betas) {
    const auto partitionMatrix = makePartitionMatrix(betas);
    VectorType result(NumBeta);
    for (size_t i = 0; i < NumBeta; ++i) {
        ScalarType sumPartition = 0;
        ScalarType numParticle = 0;
        for (unsigned int numSpinUp = 0; numSpinUp <= NumSite; ++numSpinUp) {
            for (unsigned int numSpinDown = numSpinUp; numSpinDown <= NumSite; ++numSpinDown) {
                const int factor = numSpinUp == numSpinDown ? 1 : 2;
                const ScalarType n = (numSpinUp + numSpinDown) * factor;
                sumPartition += partitionMatrix(numSpinUp, numSpinDown)[i] * factor;
                numParticle += partitionMatrix(numSpinUp, numSpinDown)[i] * n;
            }
        }
        result[i] = numParticle / (sumPartition * ScalarType(NumSite));
    }
    return result;
}

void plotTPQ() {
    const auto betas = VectorType::linspace(0, 4, NumBeta);
    const auto rhos = makeDensity(betas);

    Plot* plot = new Plot(0, 4, 0.5, 1.05, 1, 0.2);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("&beta;");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("&rho;");
    axisY->setLabelFormat("%.1f");
    plot->line(betas, rhos);
    plot->show();
    QApplication::exec();
}

void plotAll() {
    const auto betas = VectorType::linspace(0, 4, NumBeta);
    VectorType tpq_hphi, fulldiag_hphi, tpq, fulldiag;
    H5File h5f("TPQ.h5", H5File::ReadOnly);
    tpq_hphi.read(h5f, "TPQ_HPhi");
    fulldiag_hphi.read(h5f, "FullDiag_HPhi");
    tpq.read(h5f, "TPQ");
    fulldiag.read(h5f, "FullDiag");

    Plot* plot = new Plot(0, 4, 0.5, 1.05, 1, 0.2);
    plot->getLegend().setAlignment(Qt::AlignRight);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("&beta;");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("&rho;");
    axisY->setLabelFormat("%.1f");

    auto& scatter_tpq = plot->scatter(betas, tpq_hphi);
    scatter_tpq.setName("TPQ(HPhi)");
    scatter_tpq.setMarkerSize(10);
    auto pen1 = scatter_tpq.pen();
    pen1.setColor(scatter_tpq.color());
    scatter_tpq.setPen(pen1);
    scatter_tpq.setColor(Qt::transparent);

    auto& scatter_fd = plot->scatter(betas, fulldiag_hphi);
    scatter_fd.setName("FullDiag(HPhi)");
    scatter_fd.setMarkerSize(10);
    auto pen2 = scatter_fd.pen();
    pen2.setColor(scatter_fd.color());
    scatter_fd.setPen(pen2);
    scatter_fd.setColor(Qt::transparent);

    auto& line_tpq = plot->line(betas, tpq);
    line_tpq.setName("TPQ");
    line_tpq.setColor(pen1.color());

    auto& line_fd = plot->line(betas, fulldiag);
    line_fd.setName("FullDiag");
    auto pen = line_fd.pen();
    pen.setColor(pen2.color());
    pen.setStyle(Qt::DashLine);
    QList<qreal> dashes;
    dashes << 20 << 20;
    pen.setDashPattern(dashes);
    line_fd.setPen(pen);

    plot->show();
    QApplication::exec();
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    plotTPQ();
    plotAll();
    return 0;
}
