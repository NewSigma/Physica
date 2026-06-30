/*
 * Copyright 2021-2026 Weibo He.
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
#include <algorithm>
#include <iostream>
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Matrix/DenseMatrix.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/MatrixDecomp/Cholesky.h"
#include "Physica/Core/Math/Algebra/LinearAlgebra/Eigen/EigenSolver.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using T = float64;

namespace {
    class InfiniteDeepWell {
        constexpr static size_t NumBasis = 8;
        using MatrixType = DenseMatrix<T, MatrixMajor::Col, NumBasis, NumBasis>;
    public:
        static void run() {
            MatrixType overlap = getOverlapMatrix();
            //Workaround for generalised eigenvalue problem
            MatrixType cholesky = Cholesky(overlap);
            MatrixType inv_cholesky = cholesky.inv();
            MatrixType hamilton = getHamiltonMatrix();
            MatrixType hamilton_mod = (inv_cholesky * hamilton).compute() * inv_cholesky.transpose();
            EigenSolver<T> solver(hamilton_mod, true);
            const auto& eigenvalues = solver.getEigenvalues();
            std::cout << "\tNumerical\tAnalytical\n";
            Array<double, NumBasis> energy{};
            for (size_t i = 0; i < NumBasis; ++i)
                energy[i] = double(eigenvalues[i].real());
            std::ranges::sort(energy);

            for (size_t i = 0; i < NumBasis; ++i) {
                const double temp = double(i + 1) * M_PI;
                std::cout << '\t' << energy[i] << "\t\t" << (temp * temp * 0.25) << '\n';
            }

            Plot* numerical = new Plot(-1, 1, -1.2, 1.2, 0.5, 0.5);
            numerical->getChart()->setTitle("Numerical solution");
            numerical->getChart()->legend()->setAlignment(Qt::AlignRight);
            for (size_t i = 0; i < 8; ++i)
                plotWave(*numerical, solver, i);
            numerical->show();

            Plot* analytical = new Plot(-1, 1, -1.2, 1.2, 0.5, 0.5);
            analytical->getChart()->setTitle("Analytical solution");
            analytical->getChart()->legend()->setAlignment(Qt::AlignRight);
            for (size_t i = 0; i < 8; ++i)
                plotReferenceWave(*analytical, i);
            analytical->show();
            QApplication::exec();
        }
    private:
        static MatrixType getHamiltonMatrix() {
            return MatrixType::generate(NumBasis, NumBasis, [](size_t i, size_t j) -> T {
                if (j % 2 != i % 2)
                    return 0;
                const T sum = T(i + j);
                const T pro = T(i * j);
                const T numerator = T(1) - sum - T(2) * pro;
                const T denominator = (sum + T(3)) * (sum + T(1)) * (sum - T(1));
                return T(-8) * numerator / denominator;
            });
        }

        static MatrixType getOverlapMatrix() {
            return MatrixType::generate(NumBasis, NumBasis, [](size_t i, size_t j) -> T {
                if (j % 2 != i % 2)
                    return 0;
                const T sum = T(i + j);
                const T term1 = T(2) * (reciprocal(sum + T(1)) + reciprocal(sum + T(5)));
                const T term2 = T(4) * reciprocal(sum + T(3));
                return term1 - term2;
            });
        }

        static T baseFunction(const T& s, size_t n) {
            return pow(s, T(n)) * (square(s) - T(1));
        }

        static void plotWave(Plot& plot, const EigenSolver<T>& solver, size_t n) {
            constexpr size_t SampleCount = 100;
            DenseVector<T, SampleCount> x{};
            DenseVector<T, SampleCount> y{};
            const T step = T(2) / T(SampleCount);
            T temp_x = -T(1);
            for (size_t i = 0; i < SampleCount; ++i) {
                x[i] = temp_x;
                T temp_y = T(0);
                for (size_t j = 0; j < NumBasis; ++j)
                    temp_y += solver.getEigenvectors().col(n)[j].real() * baseFunction(temp_x, j);
                y[i] = temp_y;
                temp_x += step;
            }
            QString name = QString("E = %1").arg(solver.getEigenvalues()[n].real().toMachine());
            auto& spline = plot.spline(x, y);
            spline.setName(name);
        }

        static void plotReferenceWave(Plot& plot, size_t n) {
            constexpr size_t SampleCount = 100;
            DenseVector<T, SampleCount> x{};
            DenseVector<T, SampleCount> y{};
            const T step = T(2) / T(SampleCount);
            T temp = -T(1);
            const T factor = square(T(n) * T(M_PI * 0.25));
            if (n % 2U == 0) {
                for (size_t i = 0; i < SampleCount; ++i) {
                    x[i] = temp;
                    y[i] = cos(temp * factor);
                    temp += step;
                }
            }
            else {
                for (size_t i = 0; i < SampleCount; ++i) {
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
        constexpr static size_t NumBasis = 4;
        constexpr static Array<double, NumBasis> BasisCoeff{13.00773, 1.962079, 0.444529, 0.1219492};
        using MatrixType = DenseMatrix<T, MatrixMajor::Col, NumBasis, NumBasis>;
    public:
        void run() {
            MatrixType overlap = getOverlapMatrix();
            //Workaround for generalised eigenvalue problem
            MatrixType cholesky = Cholesky(overlap);
            MatrixType inv_cholesky = cholesky.inv();
            MatrixType hamilton = getHamiltonMatrix();
            MatrixType hamilton_mod = (inv_cholesky * hamilton).compute() * inv_cholesky.transpose();
            EigenSolver<T> solver(hamilton_mod, true);
            const auto& eigenvalues = solver.getEigenvalues();
            size_t groundStateIndex = eigenvalues.reals().argmin();
            for (size_t i = 0; i < NumBasis; ++i)
                std::cout << '\t' << eigenvalues[i] << '\n';

            Plot* plot = new Plot(0, 5, 0, 0.6001, 1, 0.2);

            auto eigenvectors = solver.getEigenvectors();
            VectorND<T> real_eigenvector = eigenvectors.col(groundStateIndex).reals();
            real_eigenvector = inv_cholesky.transpose() * real_eigenvector; //Safe for in-place product
            plotWave(*plot, real_eigenvector);

            plotReferenceWave(*plot);
            plot->getAxisX()->setLabelFormat("%d");
            plot->show();
            QApplication::exec();
        }
    private:
        static MatrixType getHamiltonMatrix() {
            return MatrixType::generate(NumBasis, NumBasis, [](size_t i, size_t j) {
                const T sum = T(BasisCoeff[i] + BasisCoeff[j]);
                const T pro = T(BasisCoeff[i] * BasisCoeff[j]);
                const T kinetic = T(3) * pro * sqrt(T(M_PI) / sum) * T(M_PI) / square(sum);
                const T coulomb = T(-2) * T(M_PI) / sum;
                return kinetic + coulomb;
            });
        }

        static MatrixType getOverlapMatrix() {
            return MatrixType::generate(NumBasis, NumBasis, [](size_t i, size_t j) {
                const T sum = T(BasisCoeff[i] + BasisCoeff[j]);
                const T temp = T(M_PI) / sum;
                return temp * sqrt(temp);
            });
        }

        static T baseFunction(const T& s, size_t n) {
            return exp(T(-BasisCoeff[n]) * square(s));
        }

        void plotWave(Plot& plot, const Vector auto& coeff) {
            constexpr size_t SampleCount = 100;
            DenseVector<T, SampleCount> x{};
            DenseVector<T, SampleCount> y{};
            const T step = T(5) / T(SampleCount);
            T temp_x = T(0);
            for (size_t i = 0; i < SampleCount; ++i) {
                x[i] = temp_x;
                T temp_y = T(0);
                for (size_t j = 0; j < NumBasis; ++j)
                    temp_y += coeff[j] * baseFunction(temp_x, j);
                y[i] = temp_y;
                temp_x += step;
            }
            auto& spline = plot.spline(x, y);
            spline.setName("Numerical");
        }

        static void plotReferenceWave(Plot& plot) {
            constexpr size_t SampleCount = 100;
            DenseVector<T, SampleCount> x{};
            DenseVector<T, SampleCount> y{};
            const T step = T(5) / T(SampleCount);
            T temp = T(0);
            for (size_t i = 0; i < SampleCount; ++i) {
                x[i] = temp;
                y[i] = exp(-temp) / sqrt(T(M_PI)); //The wave function in [1] is not normalized
                temp += step;
            }
            auto& spline = plot.spline(x, y);
            spline.setName("Analytical");
        }
    };

    class HeliumAtom {
        constexpr static size_t NumBasis = 4;
        constexpr static Array<double, NumBasis> BasisCoeff{0.298073, 1.242567, 5.782948, 38.474970};
    public:
        using MatrixType = DenseMatrix<T, MatrixMajor::Col, NumBasis, NumBasis>;
    public:
        void run(VectorND<T>& trial_solution, const T& criteria) {
            MatrixType overlap = getOverlapMatrix();
            //Workaround for generalised eigenvalue problem
            MatrixType cholesky = Cholesky(overlap);
            MatrixType inv_cholesky = cholesky.inv();
            auto real_eigenvalues = VectorND<T>{};
            bool stop = false;
            while (!stop) {
                MatrixType hamilton = getHamiltonMatrix(trial_solution);
                MatrixType hamilton_mod = (inv_cholesky * hamilton).compute() * inv_cholesky.transpose();
                EigenSolver<T> solver(hamilton_mod, true);

                real_eigenvalues = solver.getEigenvalues().reals();
                size_t groundStateIndex = real_eigenvalues.argmin();
                auto eigenvectors = solver.getEigenvectors();
                VectorND<T> real_eigenvector = eigenvectors.col(groundStateIndex).reals();
                real_eigenvector = inv_cholesky.transpose() * real_eigenvector; //Safe for in-place product

                stop = VectorND<T>(real_eigenvector - trial_solution).norm() < criteria;
                trial_solution = real_eigenvector;
            }
            std::cout << "Ground state energy: " << groundStateEnergy(trial_solution) << '\n';

            Plot* plot = new Plot(0, 5, 0, 1.5, 1, 0.5);
            plotWave(*plot, trial_solution);
            plot->show();
            QApplication::exec();
        }
    private:
        static T baseFunction(const T& s, size_t n) {
            return exp(T(-BasisCoeff[n]) * square(s));
        }

        static MatrixType getHamiltonMatrix(const VectorND<T>& trial_solution) {
            return MatrixType::generate(NumBasis, NumBasis, [&](size_t i, size_t j) {
                T coulomb2 = T(0);
                for (size_t m = 0; m < NumBasis; ++m)
                    for (size_t n = 0; n < NumBasis; ++n)
                        coulomb2 += calcValueQ(i, m, j, n) * trial_solution[m] * trial_solution[n];
                return calcValueH(i, j) + coulomb2;
            });
        }

        static MatrixType getOverlapMatrix() {
            return MatrixType::generate(NumBasis, NumBasis, [](size_t i, size_t j) {
                const T sum = T(BasisCoeff[i] + BasisCoeff[j]);
                const T temp = T(M_PI) / sum;
                return temp * sqrt(temp);
            });
        }

        static T calcValueH(size_t p, size_t r) {
            const T sum = T(BasisCoeff[p] + BasisCoeff[r]);
            const T pro = T(BasisCoeff[p] * BasisCoeff[r]);
            const T kinetic = T(3) * pro * sqrt(T(M_PI) / sum) * T(M_PI) / square(sum);
            const T coulomb1 = T(-4 * M_PI) / sum;
            return kinetic + coulomb1;
        }

        static T calcValueQ(size_t p, size_t r, size_t q, size_t s) {
            const T sum = T(BasisCoeff[p] + BasisCoeff[q]);
            const T sum1 = T(BasisCoeff[r] + BasisCoeff[s]);
            const T demoninator = sum * sum1 * sqrt(sum + sum1);
            const T numerator = T(2) * sqrt(T(M_PI)) * square(T(M_PI));
            return numerator / demoninator;
        }

        static T groundStateEnergy(const VectorND<T>& solution) {
            T energy = T(0);
            for (size_t i = 0; i < NumBasis; ++i) {
                for (size_t j = 0; j < NumBasis; ++j) {
                    energy += T(2) * solution[i] * solution[j] * calcValueH(i, j);
                    for (size_t m = 0; m < NumBasis; ++m)
                        for (size_t n = 0; n < NumBasis; ++n)
                            energy += calcValueQ(i, m, j, n) * solution[i] * solution[j] * solution[m] * solution[n];
                }
            }
            return energy;
        }

        void plotWave(Plot& plot, const Vector auto& coeff) {
            constexpr size_t SampleCount = 100;
            DenseVector<T, SampleCount> x{};
            DenseVector<T, SampleCount> y{};
            const T step = T(5) / T(SampleCount);
            T temp_x = T(0);
            for (size_t i = 0; i < SampleCount; ++i) {
                x[i] = temp_x;
                T temp_y = T(0);
                for (size_t j = 0; j < NumBasis; ++j)
                    temp_y += coeff[j] * baseFunction(temp_x, j);
                y[i] = temp_y;
                temp_x += step;
            }
            auto& spline = plot.spline(x, y);
            spline.setName("Numerical");
        }
    };
}
/**
 * Reference:
 * [1] J. H. Thijssen. Computational Physics[M]. London: Cambridge University Press, 2013:29-51
 */
int main(int argc, char** argv) {
    QApplication app(argc, argv);

    std::cout << "Example 1:\n";
    InfiniteDeepWell().run();

    std::cout << "Example 2:\n";
    HedrogenAtom().run();

    std::cout << "Example 3:\n";
    VectorND<T> solution{0.25, 0.25, 0.25, 0.25};
    HeliumAtom().run(solution, 1E-8);
    return 0;
}
