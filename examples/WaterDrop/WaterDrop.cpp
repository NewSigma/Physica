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
#include "Physica/Core/Math/Calculus/ODE/ODESolver.h"
#include "Physica/Core/Math/Calculus/Function/FindRoot/Bisection.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;

using T = Scalar<Double, false>;
using ODE = ODESolver<T, Vector<T>>;

struct WaterDropArgs {
    T radius;
    T sigma;
    T rho;
    T tangent;
    T g;
};

class WaterDropSolver {
private:
    ODE solver;
    T stepSize;
    T const1;
    T const2;
    T const3;
public:
    WaterDropSolver(const WaterDropArgs& drop, T rmin, T stepSize_);
    /* Operations */
    T solve();
    int output();
    /* Getters */
    T getVolume() const;
    T getMinTangent() const;
    T getLambda() const { return const3 / const2; }
private:
    void setLambda(const T& lambda) { const3 = const2 * lambda; }
};

WaterDropSolver::WaterDropSolver(const WaterDropArgs& drop, T rmin, T stepSize_)
        : solver(-drop.radius, -rmin, stepSize_, {0, drop.tangent})
        , stepSize(stepSize_)
        , const1(drop.rho * drop.g / drop.sigma)
        , const2(drop.rho / drop.sigma)
        , const3(0) {}

T WaterDropSolver::solve() {
    return bisection([&](const T& lambda) { //We use bisection method by experience and without any prove
               setLambda(lambda);
               solver.reset();
               solver.rungeKutta4([&](const T& r, const Vector<T>& v) {
                       const T momentum = v[1];
                       const T momentum_2_1 = momentum * momentum + T(1);
                       const T sqrt_momentum_2_1 = sqrt(momentum_2_1);
                       const T term1 = (momentum_2_1 * sqrt_momentum_2_1 * (const1 * v[0] + const3));
                       const T term2 = momentum_2_1 * momentum / r;
                       return Vector<T>{momentum, term1 - term2};
                   });
               return getMinTangent();
           }, T(0), T(-10), T(10)); //10 is selected by experience
}

int WaterDropSolver::output() {
    const auto& r = solver.getX();
    const auto& solution = solver.getSolution();

    const size_t length = solution.getColumn();
    Vector<T> z{};
    z.resize(length);

    T volumeHelper(0);
    for (size_t i = 0; i < length; ++i) {
        T temp1 = solution(0, i);
        z[i] = temp1;
        volumeHelper += temp1 * r[i];
    }
    std::cout << "Volume is " << abs(volumeHelper * stepSize) << " m^3" << std::endl;
    std::cout << "Minimum tangent is " << solution(1, length - 1) << std::endl;

    Plot* r_z = new Plot(r, z);
    r_z->show();
    return QApplication::exec();
}

T WaterDropSolver::getVolume() const {
    const auto& r = solver.getX();
    const auto& solution = solver.getSolution();

    const size_t length = solution.getColumn();
    T volumeHelper(0);
    for (size_t i = 0; i < length; ++i)
        volumeHelper += solution(0, i) * r[i];
    return abs(volumeHelper * stepSize);
}

T WaterDropSolver::getMinTangent() const {
    const auto& solution = solver.getSolution();
    return solution(1, solution.getColumn() - 1);
}
/**
 * This program calculates the surface equation of a water drop on a infinite smooth plain.
 * Specifying any two of mass, radius or contact angle should give the other value.
 */
int main(int argc, char** argv) {
    QApplication app(argc, argv);

    const T rmin = 0.000001;
    const T stepSize = rmin;
    const T sigma = 0.074;
    const T rho = 1000;
    const T g = 9.8;
    T radius = 0.0001;

    const size_t length = 50;
    Vector<T> r_arr{};
    r_arr.resize(length);
    Vector<T> lambda_arr{};
    lambda_arr.resize(length);
    Vector<T> volume_arr{};
    volume_arr.resize(length);
    for (size_t i = 0; i < length; ++i) {
        r_arr[i] = radius;
        WaterDropSolver solver({radius, sigma, rho, 1, 9.8}, rmin, stepSize);
        const T lambda = solver.solve();
        lambda_arr[i] = lambda;
        volume_arr[i] = solver.getVolume();
        radius += T(0.0001);
    }
    Plot* r_lambda = new Plot(r_arr, lambda_arr);
    r_lambda->show();
    Plot* r_volume = new Plot(r_arr, volume_arr);
    r_volume->show();

    WaterDropSolver solver({radius, sigma, rho, 1, 9.8}, rmin, stepSize);
    const T lambda = solver.solve();
    std::cout << "Lambda is " << lambda << std::endl;
    return solver.output();
}
