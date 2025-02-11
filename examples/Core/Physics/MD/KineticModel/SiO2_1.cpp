/*
 * Copyright 2023-2025 Weibo He. All rights reserved.
 *
 * This file is part of PhysicaNotes.
 */
#include <iostream>
#include <fstream>
#include <QApplication>
#include <QBoxLayout>
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Physica/Core/Physics/MD/ForceModel/BKSModel.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/KineticModel/FireModel.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using namespace Physica::Gui;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using KineticModel = FreeModel<ScalarType, 3, 1, RPMDIntegrator::Exact>;
using ForceModel = BKSModel<ScalarType, Ewald<ScalarType, RSpaceEwald<ScalarType, true>>, true>;
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15);
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(10);

MDCell<ScalarType> makeSystem() {
    std::ifstream fin("SiO2.vasp");
    if (!fin) {
        std::cout << "[Error]: Poscar not found" << std::endl;
        exit(EXIT_FAILURE);
    }
    Poscar<ScalarType> poscar{};
    fin >> poscar;
    poscar.scale(PhyConst<AU>::angstormToBohr(1));
    return MDCell<ScalarType>(poscar);
}

int main(int argc, char** argv) {
    RPMD<ScalarType, 3, 1> rpmd(makeSystem(), 1, 1, 0, timeStep);

    ForceModel forceModel(rpmd.phaseToCell(0), pair_cutoff, {});
    KineticModel kineticModel(0, 1);
    FireModel<ScalarType, 3> fire(timeStep, 10 * timeStep);

    VectorType f2norm(2000);
    for (size_t i = 0; i < f2norm.getLength(); ++i) {
        //rpmd.setTimeStep(fire.getTimeStep()); // Adaptive time step
        rpmd.fire_vstep<KineticModel, ForceModel, SeqExecutor>(fire, kineticModel, forceModel);
        f2norm[i] = fire.getForceNorm();
    }

    QApplication app(argc, argv);
    {
        Plot* plot = new Plot(0, 2000, -16, 0, 500, 5);
        plot->getChart()->legend()->setVisible(false);
        auto* axisX = plot->getAxisX();
        auto* axisY = plot->getAxisY();
        axisX->setTitleText("Step");
        axisX->setLabelFormat("%d");
        axisY->setTitleText("Force/a.u.");
        axisY->setLabelFormat("10<sup>%d</sup>");
        const VectorType logNorm = ln(f2norm) * reciprocal(ln(ScalarType(10)));
        plot->line(logNorm);
        plot->show();
    }
    return QApplication::exec();
}
