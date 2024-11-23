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
#include <iostream>
#include <QApplication>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Statistics/NumCharacter.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/SpinRepr.h"
#include "Physica/Core/Physics/ManyBody/TPQ.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using RandomType = Random<MT19937, 10000>;
constexpr unsigned int NumSite = 4;
constexpr double HoppingT = 1;
constexpr double RepelU = 8;
constexpr unsigned int NumSample = 100;
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
        psi.template random_normal<RandomType>(1);
        for (unsigned int j = 0; j < NumBeta; ++j) {
            if (j != 0)
                psi.nvt_step(hamilton, deltaBeta);
            toNextMean(result[j], i, psi.calcPartitionXi());
        }
    }
    return result;
}

HalfDenseMatrixStorage<VectorType> makePartitionMatrix(const VectorType& betas) {
    HalfDenseMatrixStorage<VectorType> result(NumSite + 1);
    for (unsigned int numSpinUp = 0; numSpinUp <= NumSite; ++numSpinUp) {
        for (unsigned int numSpinDown = numSpinUp; numSpinDown <= NumSite; ++numSpinDown) {
            const bool useInversionSymm = numSpinUp == numSpinDown;
            VectorType elem;
            if (useInversionSymm) {
                using ReprType = SpinRepr<1, NumSite, true>;
                elem = calcPartition(ReprType(numSpinUp, numSpinDown), betas);
            }
            else {
                using ReprType = SpinRepr<1, NumSite, false>;
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

int main(int argc, char** argv) {
    const auto betas = VectorType::linspace(0, 4, NumBeta);
    const auto rhos = makeDensity(betas);

    QApplication app(argc, argv);
    Plot* plot = new Plot(0, 4, 0.5, 1.05, 1, 0.2);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("&beta;");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("&rho;");
    axisY->setLabelFormat("%.1f");
    plot->line(betas, rhos);
    plot->show();
    return QApplication::exec();
}
