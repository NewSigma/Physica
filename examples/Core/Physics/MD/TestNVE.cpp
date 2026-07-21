/*
 * Copyright 2024-2026 Weibo He.
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
#include <QApplication>
#include <QtCharts/QLegend>
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Physica/Core/Physics/MD/ForceModel/BKSModel.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"

import Physica.Gui.Plot;

using namespace Physica;
using T = float64;
using VectorType = VectorND<T>;
using RandomSource = Random<>;
using MDType = RPMD<T, 3, 1>;
using KineticModel = FreeModel<T, 3, 1, RPMDIntegrator::Exact>;
using ForceModel = BKSModel<T, Ewald<T, RSpaceEwald<T, true>>, false>;
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15);
constexpr double pair_cutoff = PhyConst<AU>::angstromToBohr(10);
constexpr double temperatureT = PhyConst<AU>::kToTemperature(100);

namespace {
    /**
     * Initial structure from [1], modifying lattice according to [2]
     *
     * Reference:
     * [1] Materials Project, mp-7000; https://doi.org/10.17188/1272685
     * [2] Comput. Mater. Sci. 175, 109584 (2020); https://doi.org/10.1016/j.commatsci.2020.109584
     */
    MDCell<T> makeSystem() {
        using CrystalCellType = CrystalCell<T>;
        using LatticeMatrix = MDCell<T>::LatticeMatrix;
        using PositionMatrix = MDCell<T>::PositionMatrix;
        const LatticeMatrix lattice{
            5.0000000000, 0.0000000000, 0.0000000000,
        -2.5000000000, 4.3333333333, 0.0000000000,
            0.0000000000, 0.0000000000, 5.3333333333
        };
        PositionMatrix pos {
            0.256094009, 0.414853990, 0.794543028,
            0.585146010, 0.841239989, 0.127876326,
            0.158759996, 0.743906021, 0.461209685,
            0.743906021, 0.158759996, 0.538790345,
            0.414853990, 0.256094009, 0.205457002,
            0.841239989, 0.585146010, 0.872123659,
            0.523694992, 0.523694992, 0.000000000,
            0.000000000, 0.476305008, 0.666666687,
            0.476305008, 0.000000000, 0.333333343
        };

        CrystalCellType cell({lattice, pos, CrystalCellType::Type::Direct}, {8, 8, 8, 8, 8, 8, 14, 14, 14});
        cell.scale(PhyConst<AU>::angstromToBohr(1));
        MDCell<T> cell1(std::move(cell));
        cell1.toSuperCell<ExtendCellOption::AtomMajor>(4, 6, 3);
        /* To orthogonal */ {
            LatticeMatrix latt = cell1.getLattice();
            latt[1, 0] = T(0);
            cell1.setLattice(std::move(latt));
            cell1.normalize();
        }
        return cell1;
    }

    template<ExecutePolicy P>
    VectorND<T> testNVE(const auto& rpmd_, T duration, auto& kineticModel, auto& forceModel) {
        auto rpmd = rpmd_;
        const uint64_t step = rpmd.durationToStep(duration, timeStep);
        VectorND<T> ene(step);
        for (uint64_t i = 0; i < step; ++i) {
            rpmd.template nve_step<P>(kineticModel, forceModel);
            ene[i] = rpmd.calcKineticClassical() + rpmd.template calcPotential<P>(forceModel);
        }
        ene -= T(ene[0]);
        return ene;
    }
}


int main(int argc, char** argv) {
    ThreadPool::numThreadRequired = 4;
    MDType rpmd(makeSystem(), 1, 1, temperatureT, timeStep);
    rpmd.initMomentum<KineticModel, RandomSource>();
    KineticModel kineticModel(temperatureT, 1);
    ForceModel forceModel(rpmd.phaseToCell(0), pair_cutoff, {});
    const VectorType energy1 = testNVE<Thread>(rpmd, PhyConst<AU>::secondToTime(1E-12), kineticModel, forceModel);

    rpmd.setTimeStep(timeStep * 2);
    const VectorType energy2 = testNVE<Thread>(rpmd, PhyConst<AU>::secondToTime(1E-12), kineticModel, forceModel);

    rpmd.setTimeStep(timeStep * 5);
    const VectorType energy5 = testNVE<Thread>(rpmd, PhyConst<AU>::secondToTime(1E-12), kineticModel, forceModel);

    QApplication app(argc, argv);
    {
        Plot* plot = new Plot(0, 1, -0.1, 0.04, 0.2, 0.1);
        plot->getChart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
        auto* axisX = plot->getAxisX();
        auto* axisY = plot->getAxisY();
        axisX->setTitleText("Time/ps");
        axisX->setLabelFormat("%.1f");
        axisY->setTitleText("Internal energy/a.u.");
        axisY->setLabelFormat("%.1f");
        {
            const VectorType t = VectorType::linspace(0, 1, energy1.getLength());
            plot->line(t, energy1).setName("1 fs");
        }
        {
            const VectorType t = VectorType::linspace(0, 1, energy2.getLength());
            plot->line(t, energy2).setName("2 fs");
        }
        {
            const VectorType t = VectorType::linspace(0, 1, energy5.getLength());
            plot->line(t, energy5).setName("5 fs");
        }
        plot->show();
    }
    return QApplication::exec();
}
