/*
 * Copyright 2020 WeiBo He.
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
#include <memory>
#include <QtWidgets/QApplication>
#include <Physica/Core/Math/Calculus/ODE/ODESolver.h>
#include <Physica/Core/Math/Calculus/SpetialFunctions.h>
#include <Physica/Core/Physics/PhyConst.h>
#include <Physica/Gui/Plot/Plot.h>

using namespace Physica::Core;
using namespace Physica::Gui;
using T = Scalar<Double, false>;
using ODE = ODESolver<T>;

struct ProgramArgs {
    T epsilon;
    T rho;
    T reducedMass;
    T r_max;
    T energy;
    T radialNum;
    T wave_length;

    ProgramArgs(T epsilon_, T rho_, T reducedMass_, T r_max_, T energy_, T radialNum_)
            : epsilon(epsilon_)
            , rho(rho_)
            , reducedMass(reducedMass_)
            , r_max(r_max_)
            , energy(energy_)
            , radialNum(radialNum_) {
        wave_length = sqrt(T(2) / (T(reducedMass) * energy * PhyConst<SI>::unitCharge * 1E-3)) * T(M_PI * PhyConst<SI>::reducedPlanck * 1E10); //1E10 for conversion from meter to angstrom
    }
};

std::unique_ptr<ODE> solveEqu(const ProgramArgs& args) {
    const T m2_h_2 = T(2 * args.reducedMass.getTrivial() / (PhyConst<SI>::reducedPlanck * PhyConst<SI>::reducedPlanck) * PhyConst<SI>::unitCharge * 1E-23) * square(args.rho); //1E-23 for unit conversion
    const T stepSize(0.001);
    T crossSection = 0;

    const T r0 = args.rho * T(0.6);
    const T r0_2 = square(r0);
    const T r0_5 = square(r0_2) * r0;
    const T r0_6 = r0_5 * r0;
    const T c = sqrt(args.epsilon * m2_h_2 / T(25));
    ODE& solver = *(new ODE(r0 / args.rho, (args.r_max + args.wave_length / T(2)) / args.rho, stepSize, {exp(-c / r0_5)}));

    solver.degenerate_numerov([&](T h) -> T {
        const T frac = reciprocal(h);
        const T frac_2 = square(frac);
        const T frac_4 = square(frac_2);
        const T frac_6 = frac_4 * frac_2;
        const T V = args.epsilon * frac_6 * (frac_6 - T(2)); //Modified Lennard-Jones potential
        return m2_h_2 * (V - args.energy) + T(args.radialNum * (args.radialNum + T(1))) / square(h);
    }, {T(5) * c * exp(-c / r0_5) / r0_6});
    return std::unique_ptr<ODE>(&solver);
}

T calcPhase(const ProgramArgs& args, const Vector<T>& h, const Vector<T>& wave) {
    const T wave_vector = T(2 * M_PI) / args.wave_length;
    const size_t length = h.getLength();
    assert(length == wave.getLength());
    T r1, r2, wave1, wave2;
    /* Read values */ {
        const T h_max = args.r_max / args.rho;
        for (size_t i = 0; i < length; ++i) {
            if (h[i] > h_max) {
                r1 = h[i] * args.rho;
                wave1 = wave[i];
                break;
            }
        }
        r2 = h[length - 1] * args.rho;
        wave2 = wave[length - 1];
    }
    const T factor = r1 * wave2 / (r2 * wave1);
    const T k_r1 = wave_vector * r1;
    const T k_r2 = wave_vector * r2;
    T j1, j2, n1, n2, unused1, unused2;
    sphericalBesselJn_Yn_dJn_dYn(args.radialNum, k_r1, j1, n1, unused1, unused2);
    sphericalBesselJn_Yn_dJn_dYn(args.radialNum, k_r2, j2, n2, unused1, unused2);
    return arctan(T(factor * j1 - j2) / T(factor * n1 - n2));
}

T calcCrossSection(double energy) {
    const double reducedMass = 83.8 / (83.8 + 1) * PhyConst<SI>::protonMass; //83.8 is the relative atomic mass of Kr
    ProgramArgs args(5.9, 3.57, reducedMass, 5 * 3.57, energy, 0);
    const T l_max = T(2 * M_PI) * args.r_max / args.wave_length;
    T crossSection = 0;

    for (; args.radialNum < l_max; args.radialNum += T(1)) {
        auto solver = solveEqu(args);
        const auto& h = solver->getX();
        const auto& solution = solver->getSolution();

        Vector<T> waveFunc(solution.getColumn());
        for (size_t i = 0; i < solution.getColumn(); ++i)
            waveFunc[i] = solution[i][0];
        const T phase = calcPhase(args, h, waveFunc);
        crossSection += (T(2) * args.radialNum + T(1)) * square(sin(phase));
    }
    return crossSection * square(args.wave_length);
}

void plotPWBaseWave(double energy, double radialNum) {
    const double reducedMass = 83.8 / (83.8 + 1) * PhyConst<SI>::protonMass; //83.8 is the relative atomic mass of Kr
    ProgramArgs args(5.9, 3.57, reducedMass, 5 * 3.57, energy, radialNum);
    auto solver = solveEqu(args);
    const auto& h = solver->getX();
    const auto& solution = solver->getSolution();

    Vector<T> waveFunc(solution.getColumn());
    for (size_t i = 0; i < solution.getColumn(); ++i)
        waveFunc[i] = solution[i][0] / (h[i] * args.rho);
    Plot* plot = new Plot();
    plot->spline(h, waveFunc);
    auto& chart = *plot->chart();
    chart.setTitle("h-wave function");
    chart.legend()->hide();
    auto axes = chart.axes();
    chart.createDefaultAxes();
    chart.axes(Qt::Horizontal).first()->setTitleText("h");
    chart.axes(Qt::Vertical).first()->setTitleText("wave function");
    plot->show();
}
/**
 * Reference:
 * [1] Jos Thijssen. Computational Physics[M].London: Cambridge university press, 2013:14-28
 * [2] J.P. Toennies, W. Welz, and G. Wolf,'Molecular beam scattering studies of orbiting resonances
 * and the determination of Van der Waals potentials for H-He, Ar, Kr, and Xe and for H2-Ar, Kr and Xe,'
 * J.Chem.Phys., 71(1979), 614-42
 */
int main(int argc, char** argv) {
    const double step = 0.005;
    const double from = 0.2;
    const double to = 3.5;
    const size_t count = static_cast<size_t>((to - from) / step);
    Vector<T> energyArr(count);
    Vector<T> crossSectionArr(count);
    /* Get energy-crossSection */ {
        double temp = 0.2;
        for (size_t i = 1; i <= count; ++i) {
            energyArr[i - 1] = temp;
            crossSectionArr[i - 1] = calcCrossSection(temp) / T(3.57 * 3.57);
            temp += step;
        }
    }
    QApplication app(argc, argv);

    Plot* plot = new Plot();
    plot->spline(energyArr, crossSectionArr);
    auto& chart = *plot->chart();
    chart.setTitle("E-cross section");
    chart.legend()->hide();
    chart.createDefaultAxes();
    chart.axes(Qt::Horizontal).first()->setTitleText("E/meV");
    chart.axes(Qt::Vertical).first()->setTitleText("cross section/rho^2");
    plot->show();

    plotPWBaseWave(0.4764, 4);
    return QApplication::exec();
}