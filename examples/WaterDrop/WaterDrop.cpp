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

struct EquationEnv {
    Vector<T>& r_arr;
    Vector<T>& lambda_arr;
    Vector<T>& volume_arr;
    size_t from;
    size_t to;
    double rmin;
    double solveStepSize;
    double plotStepSize;
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
           }, T(0), T(-1.1), T(1.1)); //1.1 is selected by experience
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
    /* Plot z */ {
        Plot* r_z = new Plot();
        r_z->spline(r, z);
        auto& r_z_chart = *r_z->chart();
        r_z_chart.setTitle("r-z");
        r_z_chart.legend()->hide();
        auto r_z_axes = r_z_chart.axes();
        r_z_axes[0]->setTitleText("r/m");
        r_z_axes[1]->setTitleText("z/m");
        r_z->show();
    }
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

void parallelSolve(WaterDropArgs args, const EquationEnv& env) {
    assert(env.from < env.to);
    assert(env.to <= env.r_arr.getLength());
    assert(env.r_arr.getLength() <= env.lambda_arr.getLength());
    assert(env.r_arr.getLength() <= env.volume_arr.getLength());
    for (size_t i = env.from; i < env.to; ++i) {
        env.r_arr[i] = args.radius;
        WaterDropSolver solver(args, env.rmin, env.solveStepSize);
        const T lambda = solver.solve();
        env.lambda_arr[i] = lambda;
        env.volume_arr[i] = solver.getVolume();
        args.radius += T(env.plotStepSize);
    }
}
/**
 * This program calculates the surface equation of a water drop on a infinite smooth plain.
 * Specifying any two of mass, radius or contact angle should give the other value.
 */
int main(int argc, char** argv) {
    QApplication app(argc, argv);

    const T sigma = 0.074;
    const T rho = 1000;
    const T g = 9.8;
    const T startRadius = 0.0001;
    const double rmin = 0.000001;
    const double solveStepSize = rmin;
    const double plotStepSize = 0.0001;

    const size_t length = 50;
    Vector<T> r_arr{};
    r_arr.resize(length);
    Vector<T> lambda_arr{};
    lambda_arr.resize(length);
    Vector<T> volume_arr{};
    volume_arr.resize(length);
    {
        const unsigned int threadCount = 5;
        const unsigned int taskPerThread = length / threadCount;
        const unsigned int halfTaskPerThread = taskPerThread / 2U;
        assert(length % threadCount == 0);
        assert(taskPerThread % 2 == 0);
        std::thread threads[5];
        for (unsigned int i = 0; i < threadCount; ++i) {
            threads[i] = std::thread([&, i = i]() {
                WaterDropArgs args{startRadius, sigma, rho, 1, 9.8};
                args.radius = startRadius + T(halfTaskPerThread * i * plotStepSize);
                parallelSolve(args, { r_arr
                                    , lambda_arr
                                    , volume_arr
                                    , i * halfTaskPerThread
                                    , (i + 1) * halfTaskPerThread
                                    , rmin
                                    , solveStepSize
                                    , plotStepSize });
                args.radius = startRadius + T((length - (i + 1) * halfTaskPerThread) * plotStepSize);
                parallelSolve(args, { r_arr
                                    , lambda_arr
                                    , volume_arr
                                    , length - (i + 1) * halfTaskPerThread
                                    , length - i * halfTaskPerThread
                                    , rmin
                                    , solveStepSize
                                    , plotStepSize });
            });
        }
        for (unsigned int i = 0; i < threadCount; ++i)
            threads[i].join();
    }
    /* Plot lambda */ {
        Plot* r_lambda = new Plot();
        r_lambda->spline(r_arr, lambda_arr);
        auto& r_lambda_chart = *r_lambda->chart();
        r_lambda_chart.legend()->hide();
        r_lambda_chart.setTitle("r-lambda");
        auto r_lambda_axes = r_lambda_chart.axes();
        r_lambda_axes[0]->setTitleText("r/m");
        r_lambda_axes[1]->setTitleText("lambda");
        r_lambda->show();
    }
    /* Plot volume */ {
        Plot* r_volume = new Plot();
        r_volume->spline(r_arr, volume_arr);
        auto& r_volume_chart = *r_volume->chart();
        r_volume_chart.legend()->hide();
        r_volume_chart.setTitle("r-volume");
        auto r_volume_axes = r_volume_chart.axes();
        r_volume_axes[0]->setTitleText("r/m");
        r_volume_axes[1]->setTitleText("volume/m^3");
        r_volume->show();
    }

    WaterDropSolver solver({startRadius, sigma, rho, 1, 9.8}, rmin, solveStepSize);
    const T lambda = solver.solve();
    std::cout << "Lambda is " << lambda << std::endl;
    return solver.output();
}
