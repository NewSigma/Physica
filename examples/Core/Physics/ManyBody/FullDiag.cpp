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
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/FermiRepr.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using ScalarType = float64; // Will overflow if use float32
using VectorType = VectorND<ScalarType>;
using RandomType = Random<MT19937>;
using VectorPair = std::pair<VectorType, VectorType>;
using SpectrumMatrix = HalfDenseMatrixStorage<VectorPair, Dynamic>;
constexpr unsigned int NumSite = 4;
constexpr double HoppingT = 1;
constexpr double RepelU = 8;

template<class ReprType>
VectorPair calcPair(ReprType repr_) {
    using Hamilton = HubbardMatrix<ScalarType, ReprType>;
    const LatticeModel<1> lattice({NumSite}, 1);
    const Hubbard<ScalarType, 1> hubbard(lattice, HoppingT, RepelU);
    const Hamilton hamilton(hubbard, std::move(repr_));
    const size_t numState = hamilton.getNumState();
    SymmEigenSolver<ScalarType> solver(numState);
    solver.compute(DenseSymmMatrix<ScalarType>(hamilton), true);

    const auto eigenvectors = solver.getEigenvectors();
    const auto& repr = hamilton.getRepr();
    VectorType numParticles(numState);
    for (size_t i = 0; i < numState; ++i)
        numParticles += eigenvectors.row(i).squaredNorms() * ScalarType(repr[i].getNumParticle());
    return std::make_pair(VectorType(solver.getEigenvalues()), std::move(numParticles));
}

SpectrumMatrix makeSpectrum() {
    SpectrumMatrix result(NumSite + 1, NumSite + 1, std::make_pair(VectorType{}, VectorType{}));
    for (unsigned int numSpinUp = 0; numSpinUp <= NumSite; ++numSpinUp) {
        for (unsigned int numSpinDown = numSpinUp; numSpinDown <= NumSite; ++numSpinDown) {
            const bool isVacuum = numSpinDown == 0;
            if (isVacuum) {
                result(numSpinUp, numSpinDown) = std::make_pair<VectorType, VectorType>({0}, {0});
                continue;
            }

            const bool useInversionSymm = numSpinUp == numSpinDown;
            VectorPair pair;
            if (useInversionSymm) {
                using ReprType = FermiRepr<1, NumSite, true>;
                pair = calcPair(ReprType(numSpinUp, numSpinDown));
            }
            else {
                using ReprType = FermiRepr<1, NumSite, false>;
                pair = calcPair(ReprType(numSpinUp, numSpinDown));
            }
            result(numSpinUp, numSpinDown) = std::move(pair);
        }
    }
    return result;
}

VectorType makeDensity(const VectorType& betas) {
    const auto spectrum = makeSpectrum();
    VectorType result(betas.getLength());
    for (size_t i = 0; i < betas.getLength(); ++i) {
        const ScalarType beta = betas[i];
        ScalarType numParticles = 0, partition = 0;
        for (unsigned int numSpinUp = 0; numSpinUp <= NumSite; ++numSpinUp) {
            for (unsigned int numSpinDown = numSpinUp; numSpinDown <= NumSite; ++numSpinDown) {
                const VectorPair& pair = spectrum(numSpinUp, numSpinDown);
                VectorType weights = exp(pair.first * (-beta));
                if (numSpinDown != numSpinUp)
                    weights *= ScalarType(2);
                numParticles += pair.second * weights;
                partition += weights.sum();
            }
        }
        result[i] = numParticles / partition;
    }
    result *= reciprocal(ScalarType(NumSite));
    return result;
}

int main(int argc, char** argv) {
    const auto betas = VectorType::linspace(0, 4, 41);
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
