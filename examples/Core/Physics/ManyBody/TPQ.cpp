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
#include <Physica/Core/Math/Random/RandomPool.h>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h>
#include <Physica/Core/Math/Statistics/NumCharacter.h>
#include <Physica/Core/Physics/ManyBody/ExactDiag/Hamilton/HubbardMatrix.h>
#include <Physica/Core/Physics/ManyBody/ExactDiag/ReprSpace/SpinRepr.h>
#include <Physica/Core/Physics/ManyBody/ExactDiag/TPQ.h>
#include <Physica/Gui/Plot/Plot.h>

using namespace Physica::Core;
using namespace Physica::Gui;
using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;
using RandomPoolType = RandomPool<std::mt19937>;
constexpr unsigned int NumSite = 4;
constexpr unsigned int NumSpinUp = NumSite / 2;
constexpr unsigned int NumSpinDown = NumSite / 2;
constexpr double HoppingT = 1.0;
constexpr double repelU = 8.0;
constexpr unsigned int NumSample = 100;
constexpr unsigned int NumSystem = 6;
constexpr unsigned int NumBeta = 40;
using ReprType = SpinRepr<1, NumSite, NumSpinUp == NumSpinDown>;
using Hamilton = HubbardMatrix<ScalarType, ReprType>;

std::pair<VectorType, ScalarType> energyFullDiag(const Hamilton& hamilton, const VectorType& betas) {
    SymmEigenSolver<ScalarType> solver(hamilton.getNumState());
    solver.compute(hamilton, false);
    solver.sort();
    const auto& eigenvalues = solver.getEigenvalues();

    VectorType energys(NumBeta);
    for (size_t i = 0; i < NumBeta; ++i) {
        const VectorType weights = exp(eigenvalues * (-betas[i]));
        energys[i] = eigenvalues * weights / weights.sum();
    }
    return std::make_pair(energys, eigenvalues[0]);
}

VectorType energyTPQ(const Hamilton& hamilton, ScalarType deltaBeta) {
    auto& gen = RandomPoolType::getInstance().getGen();
    VectorType energys(NumBeta);
    auto psi = TPQ<ScalarType>::random_normal(hamilton.getNumState(), gen);
    psi.pre_nvt_step(hamilton, deltaBeta);
    for (size_t i = 0; i < NumBeta; ++i) {
        if (i != 0)
            psi.nvt_step(hamilton, deltaBeta);
        energys[i] = (psi * VectorType(hamilton * psi)) / psi.squaredNorm();
    }
    return energys;
}

int main(int argc, char** argv) {
    const LatticeModel<1> lattice({NumSite}, 1);
    const Hubbard<ScalarType, 1> hubbard(lattice, HoppingT, repelU);
    const Hamilton hamilton(hubbard, ReprType(NumSpinUp, NumSpinDown));
    const auto betas = VectorType::linspace(0, 4, NumBeta);
    const auto pair = energyFullDiag(hamilton, betas);
    const auto& energy0 = pair.first;
    const auto groundE = pair.second;

    VectorType means(NumBeta), vars(NumBeta);
    for (unsigned int sys = 0; sys < NumSystem; ++sys) {
        VectorType energys(NumBeta);
        for (unsigned int sample = 0; sample < NumSample; ++sample)
            toNextMean(energys, sample, energyTPQ(hamilton, betas[1]));
        toNextVariance(vars, means, sys, energys);
    }

    QApplication app(argc, argv);
    Plot* plot = new Plot(0, 4, -1.8, 9, 1, 2);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("&beta;");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("E");
    axisY->setLabelFormat("%d");
    {
        const VectorType devia = sqrt(vars) * ScalarType(2);
        auto& area = plot->area_center(betas, means, devia);
        const auto original = area.color();
        auto areaColor = original;
        areaColor.setAlpha(75);
        area.setColor(areaColor);

        auto& line = plot->line(betas, means);
        line.setName("TPQ");
        line.setColor(original);
    }
    {
        auto& scatter = plot->scatter(betas, energy0);
        scatter.setName("FullDiag");
        scatter.setColor(Qt::red);
        scatter.setMarkerSize(10);
    }
    {
        auto& line = plot->line(VectorType{0, 4}, VectorType{groundE, groundE});
        auto pen = line.pen();
        pen.setStyle(Qt::DashLine);
        pen.setColor(Qt::blue);
        line.setPen(pen);
    }
    plot->show();
    return QApplication::exec();
}
