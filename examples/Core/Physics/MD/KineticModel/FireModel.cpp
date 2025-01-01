/*
 * Copyright 2023-2024 Weibo He. All rights reserved.
 *
 * This file is part of PhysicaNotes.
 */
#include <iostream>
#include <QApplication>
#include "Physica/Core/Physics/MD/ForceModel/Ewald/Ewald.h"
#include "Physica/Core/Physics/MD/ForceModel/Q_TIP4P.h"
#include "Physica/Core/Physics/MD/RPMD.h"
#include "Physica/Core/Physics/MD/KineticModel/FreeModel.h"
#include "Physica/Core/Physics/MD/KineticModel/FireModel.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using namespace Physica::Gui;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;
using RandomType = Random<MT19937, 10002>;
using KineticModel = FreeModel<ScalarType, 3, 1, RPMDIntegrator::Exact>;
using ForceModel = Q_TIP4P<ScalarType, Ewald<ScalarType, RSpaceEwald<ScalarType, true>>>;
constexpr double timeStep = PhyConst<AU>::secondToTime(1E-15) * 0.5;
constexpr double pair_cutoff = PhyConst<AU>::angstormToBohr(9);

MDCell<ScalarType> makeSystem() {
    using CrystalCellType = CrystalCell<ScalarType>;
    using LatticeMatrix = MDCell<ScalarType>::LatticeMatrix;
    using PositionMatrix = MDCell<ScalarType>::PositionMatrix;
    const LatticeMatrix lattice{
         4.6635062604325164,    0.2499522611778955,    0.0000000000000000,
         2.1629745970109657,    4.1943944839773311,    0.0000000000000000,
         0.0000000000000000,    0.0000000000000000,   18.0000000000000000
    };
    PositionMatrix pos {
        0.4553508091084409, 0.3980437584135783, 0.6240303800896787,
        0.4937103263031835, 0.4030549988960055, 0.4488679230950712,
        0.5596918259357793, 0.8517822319914985, 0.6226285591691945,
        0.3686253245184842, 0.6403194088783717, 0.5163388989450929,
        0.5980496296529945, 0.8452761470000689, 0.4474667297556534,
        0.1076134753395919, 0.7454143691003363, 0.6235631952109221,
        0.6847635746074375, 0.9781822048654215, 0.5551553884196571,
        0.9457728898525665, 0.8447981726837335, 0.4479330783198919,
        0.6786065180112657, 0.9738018598142685, 0.6107906720462844,
        0.3747794285285615, 0.6414996922536187, 0.4607035572233092,
        0.7021261874138659, 0.9803177844871507, 0.4573706719037773,
        0.3512600170478342, 0.6392493670479714, 0.6141244566914832
    };
    
    CrystalCellType cell({lattice, pos, CrystalCellType::Type::Direct}, {1, 1, 1, 1, 1, 1, 1, 1, 8, 8, 8, 8});
    cell.scale(PhyConst<AU>::angstormToBohr(1));
    MDCell<ScalarType> cell1(std::move(cell));
    ForceModel::sortPosition(cell1);

    pos = cell1.getPos();
    std::normal_distribution<double> dist(0, 0.5);
    for (auto& elem : pos)
        elem += ScalarType::random_any<RandomType>(dist); //Perturbation
    cell1.setPos(pos);
    return cell1;
}

int main(int argc, char** argv) {
    RPMD<ScalarType, 3, 1> rpmd(makeSystem(), 1, 1, 0, timeStep);

    ForceModel forceModel(rpmd.phaseToCell(0), pair_cutoff, {});
    KineticModel kineticModel(0, 1);
    FireModel<ScalarType, 3> fire(timeStep, 10 * timeStep);

    VectorType f2norm(2000);
    for (size_t i = 0; i < f2norm.getLength(); ++i) {
        rpmd.fire_vstep<KineticModel, ForceModel, SequentialExecutor>(fire, kineticModel, forceModel);
        f2norm[i] = fire.getForceNorm();
    }

    const VectorType logNorm = ln(f2norm) * reciprocal(ln(ScalarType(10)));

    QApplication app(argc, argv);
    Plot* plot = new Plot(0, 2000, -16, 0, 500, 5);
    plot->getChart()->legend()->setVisible(false);
    auto* axisX = plot->getAxisX();
    auto* axisY = plot->getAxisY();
    axisX->setTitleText("Step");
    axisX->setLabelFormat("%d");
    axisY->setTitleText("Force/a.u.");
    axisY->setLabelFormat("10<sup>%d</sup>");
    plot->line(logNorm);
    plot->show();
    return QApplication::exec();
}
