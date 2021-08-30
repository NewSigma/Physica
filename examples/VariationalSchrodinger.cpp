/*
 * Copyright 2021 WeiBo He.
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
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/Transpose.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/MatrixDecomposition/Cholesky.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/EigenSolver.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;
using namespace Physica::Utils;
using ScalarType = Scalar<Double, false>;

class InfiniteDeepWell {
    constexpr static size_t baseSetCount = 8;
    using MatrixType = DenseMatrix<ScalarType,
                                DenseMatrixOption::Column | DenseMatrixOption::Vector,
                                baseSetCount,
                                baseSetCount>;
public:
    int execute(int argc, char** argv) {
        MatrixType overlap = getOverlapMatrix();
        //Workaround for generalised eigenvalue problem
        MatrixType cholesky = Cholesky(overlap);
        MatrixType inv_cholesky = cholesky.inverse();
        MatrixType hamilton = getHamiltonMatrix();
        MatrixType hamilton_mod = (inv_cholesky * hamilton).compute() * inv_cholesky.transpose();
        EigenSolver solver(hamilton_mod, true);
        const auto& eigenvalues = solver.getEigenvalues();
        std::cout << "\tNumerical\tAnalytical\n";
        std::array<double, baseSetCount> energy{};
        for (size_t i = 0; i < baseSetCount; ++i)
            energy[i] = double(eigenvalues[i].getReal());
        std::sort(energy.begin(), energy.end());

        for (size_t i = 0; i < baseSetCount; ++i) {
            const double temp = (i + 1) * M_PI;
            std::cout << '\t' << energy[i] << "\t\t" << (temp * temp * 0.25) << '\n';
        }

        QApplication app(argc, argv);
        Plot* numerical = new Plot();
        numerical->chart()->setTitle("Numerical solution");
        numerical->chart()->legend()->setAlignment(Qt::AlignRight);
        for (size_t i = 0; i < 8; ++i)
            plotWave(*numerical, solver, i);
        numerical->show();

        Plot* analytical = new Plot();
        analytical->chart()->setTitle("Analytical solution");
        analytical->chart()->legend()->setAlignment(Qt::AlignRight);
        for (size_t i = 0; i < 8; ++i)
            plotReferenceWave(*analytical, i);
        analytical->show();
        return QApplication::exec();
    }
private:
    MatrixType getHamiltonMatrix() {
        MatrixType result = MatrixType::Zeros(baseSetCount, baseSetCount);
        for (size_t i = 0; i < baseSetCount; ++i) {
            for (size_t j = i % 2; j < baseSetCount; j += 2) {
                const ScalarType sum = ScalarType(i + j);
                const ScalarType pro = ScalarType(i * j);
                const ScalarType numerator = ScalarType::One() - sum - ScalarType::Two() * pro;
                const ScalarType denominator = (sum + ScalarType(3)) * (sum + ScalarType::One()) * (sum - ScalarType::One());
                result(i, j) = ScalarType(-8) * numerator / denominator;
            }
        }
        return result;
    }

    MatrixType getOverlapMatrix() {
        MatrixType result = MatrixType::Zeros(baseSetCount, baseSetCount);
        for (size_t i = 0; i < baseSetCount; ++i) {
            for (size_t j = i % 2; j < baseSetCount; j += 2) {
                const ScalarType sum = ScalarType(i + j);
                const ScalarType term1 = ScalarType::Two() * (reciprocal(sum + ScalarType::One()) + reciprocal(sum + ScalarType(5)));
                const ScalarType term2 = ScalarType(4) * reciprocal(sum + ScalarType(3));
                result(i, j) = term1 - term2;
            }
        }
        return result;
    }

    ScalarType baseFunction(const ScalarType& s, size_t n) {
        return pow(s, ScalarType(n)) * (square(s) - ScalarType::One());
    }

    template<class MatrixType>
    void plotWave(Plot& plot, const EigenSolver<MatrixType>& solver, size_t n) {
        constexpr size_t sampleCount = 100;
        Array<ScalarType, sampleCount> x{};
        Array<ScalarType, sampleCount> y{};
        const ScalarType step = ScalarType::Two() / ScalarType(sampleCount);
        ScalarType temp_x = -ScalarType::One();
        for (size_t i = 0; i < sampleCount; ++i) {
            x[i] = temp_x;
            ScalarType temp_y = ScalarType::Zero();
            for (size_t j = 0; j < baseSetCount; ++j)
                temp_y += solver.getEigenvectors().col(n)[j].getReal() * baseFunction(temp_x, j);
            y[i] = temp_y;
            temp_x += step;
        }
        QString name = QString("E = %1").arg(solver.getEigenvalues()[n].getReal().getTrivial());
        auto& spline = plot.spline(x, y);
        spline.setName(name);
    }

    void plotReferenceWave(Plot& plot, size_t n) {
        constexpr size_t sampleCount = 100;
        Array<ScalarType, sampleCount> x{};
        Array<ScalarType, sampleCount> y{};
        const ScalarType step = ScalarType::Two() / ScalarType(sampleCount);
        ScalarType temp = -ScalarType::One();
        const ScalarType factor = square(ScalarType(n * M_PI * 0.25));
        if (n % 2U == 0) {
            for (size_t i = 0; i < sampleCount; ++i) {
                x[i] = temp;
                y[i] = cos(temp * factor);
                temp += step;
            }
        }
        else {
            for (size_t i = 0; i < sampleCount; ++i) {
                x[i] = temp;
                y[i] = sin(temp * factor);
                temp += step;
            }
        }
        auto& spline = plot.spline(x, y);
        QString name = QString("N = %1").arg(n);
        spline.setName(name);
    }
};

class HedrogenAtom {
    constexpr static size_t baseSetCount = 4;
    constexpr static double baseSetCoeff[baseSetCount]{13.00773, 1.962079, 0.444529, 0.1219492};
    using MatrixType = DenseMatrix<ScalarType,
                                DenseMatrixOption::Column | DenseMatrixOption::Vector,
                                baseSetCount,
                                baseSetCount>;
public:
    int execute(int argc, char** argv) {
        MatrixType overlap = getOverlapMatrix();
        //Workaround for generalised eigenvalue problem
        MatrixType cholesky = Cholesky(overlap);
        MatrixType inv_cholesky = cholesky.inverse();
        MatrixType hamilton = getHamiltonMatrix();
        MatrixType hamilton_mod = (inv_cholesky * hamilton).compute() * inv_cholesky.transpose();
        EigenSolver solver(hamilton_mod, true);
        const auto& eigenvalues = solver.getEigenvalues();
        size_t groundStateIndex = 0;
        for (size_t i = 0; i < baseSetCount; ++i) {
            std::cout << '\t' << eigenvalues[i] << '\n';
            if (eigenvalues[i].getReal() < eigenvalues[groundStateIndex].getReal())
                groundStateIndex = i;
        }
        QApplication app(argc, argv);
        Plot* plot = new Plot();
        plotWave(*plot, solver.getEigenvectors().col(groundStateIndex));
        plotReferenceWave(*plot);
        plot->show();
        return QApplication::exec();
    }
private:
    MatrixType getHamiltonMatrix() {
        MatrixType result = MatrixType::Zeros(baseSetCount, baseSetCount);
        for (size_t i = 0; i < baseSetCount; ++i) {
            for (size_t j = 0; j < baseSetCount; ++j) {
                const ScalarType sum = ScalarType(baseSetCoeff[i] + baseSetCoeff[j]);
                const ScalarType pro = ScalarType(baseSetCoeff[i] * baseSetCoeff[j]);
                const ScalarType kinetic = ScalarType(3) * pro * sqrt(ScalarType(M_PI) / sum) * ScalarType(M_PI) / square(sum);
                const ScalarType coulomb = ScalarType(-2) * ScalarType(M_PI) / sum;
                result(i, j) = kinetic + coulomb;
            }
        }
        return result;
    }

    MatrixType getOverlapMatrix() {
        MatrixType result = MatrixType::Zeros(baseSetCount, baseSetCount);
        for (size_t i = 0; i < baseSetCount; ++i) {
            for (size_t j = 0; j < baseSetCount; ++j) {
                const ScalarType sum = ScalarType(baseSetCoeff[i] + baseSetCoeff[j]);
                const ScalarType temp = ScalarType(M_PI) / sum;
                result(i, j) = temp * sqrt(temp);
            }
        }
        return result;
    }

    ScalarType baseFunction(const ScalarType& s, size_t n) {
        return exp(ScalarType(-baseSetCoeff[n]) * square(s));
    }

    template<class VectorType>
    void plotWave(Plot& plot, const LValueVector<VectorType>& coeff) {
        constexpr size_t sampleCount = 100;
        Array<ScalarType, sampleCount> x{};
        Array<ScalarType, sampleCount> y{};
        const ScalarType step = ScalarType::One() / ScalarType(sampleCount);
        ScalarType temp_x = ScalarType::Zero();
        for (size_t i = 0; i < sampleCount; ++i) {
            x[i] = temp_x;
            ScalarType temp_y = ScalarType::Zero();
            for (size_t j = 0; j < baseSetCount; ++j)
                temp_y += coeff[j].getReal() * baseFunction(temp_x, j);
            y[i] = temp_y;
            temp_x += step;
        }
        auto& spline = plot.spline(x, y);
        spline.setName("Numerical");
    }

    void plotReferenceWave(Plot& plot) {
        constexpr size_t sampleCount = 100;
        Array<ScalarType, sampleCount> x{};
        Array<ScalarType, sampleCount> y{};
        const ScalarType step = ScalarType::One() / ScalarType(sampleCount);
        ScalarType temp = ScalarType::Zero();
        for (size_t i = 0; i < sampleCount; ++i) {
            x[i] = temp;
            y[i] = ScalarType::Two() * exp(-temp);
            temp += step;
        }
        auto& spline = plot.spline(x, y);
        spline.setName("Analytical");
    }
};
/**
 * Reference:
 * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:29-37
 */
int main(int argc, char** argv) {
    int exit_code = 0;
    std::cout << "Example 1:\n";
    exit_code |= InfiniteDeepWell().execute(argc, argv);
    std::cout << "Example 2:\n";
    exit_code |= HedrogenAtom().execute(argc, argv);
    return exit_code;
}
