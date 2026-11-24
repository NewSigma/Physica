/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <QApplication>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseSymmMatrix.h"
#include "Physica/Core/Math/Random/Random.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/SymmEigenSolver.h"
#include "Physica/Core/Physics/ManyBody/Hamilton/HubbardMatrix.h"
#include "Physica/Core/Physics/ManyBody/ReprSpace/FermiRepr.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<>;
using VectorPair = std::pair<VectorND<T>, VectorND<T>>;
using SpectrumMatrix = SymmArray<VectorPair, Dynamic>;
constexpr int NumSite = 4;
constexpr double HoppingT = 1;
constexpr double RepelU = 8;
constexpr int MaxBeta = 16;
constexpr int NumBeta = MaxBeta * 10 + 1;

namespace {
    template<class ReprType>
    VectorPair calcPair(ReprType repr_) {
        using Hamilton = HubbardMatrix<T, ReprType>;
        const SquareLattice<1> lattice({NumSite}, 1);
        const Hamilton hamilton(HoppingT, RepelU, lattice, std::move(repr_));
        const size_t numState = hamilton.getNumState();
        SymmEigenSolver<T> solver(numState, true);
        solver.compute(DenseSymmMatrix<T>(hamilton));

        const auto eigenvectors = solver.getEigenvectors();
        const auto& repr = hamilton.getRepr();
        VectorND<T> numParticles(numState, 0);
        for (size_t i = 0; i < numState; ++i)
            numParticles += eigenvectors.col(i).squaredNorms() * T(repr[i].getNumParticle());
        return std::make_pair(VectorND<T>(solver.getEigenvalues()), std::move(numParticles));
    }

    SpectrumMatrix makeSpectrum() {
        SpectrumMatrix result(NumSite + 1, std::make_pair(VectorND<T>{}, VectorND<T>{}));
        for (int numSpinUp = 0; numSpinUp <= NumSite; ++numSpinUp) {
            for (int numSpinDown = numSpinUp; numSpinDown <= NumSite; ++numSpinDown) {
                const bool isVacuum = numSpinDown == 0;
                if (isVacuum) {
                    result[numSpinUp, numSpinDown] = std::make_pair<VectorND<T>, VectorND<T>>({0}, {0});
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
                result[numSpinUp, numSpinDown] = std::move(pair);
            }
        }
        return result;
    }

    VectorND<T> makeDensity(const VectorND<T>& betas) {
        const auto spectrum = makeSpectrum();
        VectorND<T> result(betas.getLength());
        for (size_t i = 0; i < betas.getLength(); ++i) {
            const T beta = betas[i];
            T numParticles = 0, partition = 0;
            for (int numSpinUp = 0; numSpinUp <= NumSite; ++numSpinUp) {
                for (int numSpinDown = numSpinUp; numSpinDown <= NumSite; ++numSpinDown) {
                    const VectorPair& pair = spectrum[numSpinUp, numSpinDown];
                    VectorND<T> weights = exp(pair.first * (-beta));
                    if (numSpinDown != numSpinUp)
                        weights *= T(2);
                    numParticles += pair.second * weights;
                    partition += weights.sum();
                }
            }
            result[i] = numParticles / partition;
        }
        result *= reciprocal(T(NumSite));
        return result;
    }
}

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    Plot* plot = new Plot(0, MaxBeta, 0.5, 1.05, 8, 0.2);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("&beta;");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("&rho;");
    axisY->setLabelFormat("%.1f");
    {
        const auto betas = VectorND<T>::linspace(0, MaxBeta, NumBeta);
        const auto rhos = makeDensity(betas);
        plot->line(betas, rhos);
    }
    plot->show();
    return QApplication::exec();
}
