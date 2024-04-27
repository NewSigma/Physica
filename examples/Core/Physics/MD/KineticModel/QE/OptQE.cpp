/*
 * Copyright 2023 WeiBo He. All rights reserved.
 *
 * This file is part of PhysicaNotes.
 */
#include <iostream>
#include <QApplication>
#include "Physica/Core/IO/QE/PWscfOut.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;
using namespace Physica::Gui;

using ScalarType = Scalar<Double>;
using VectorType = Vector<ScalarType>;

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    {
        Plot* plot = new Plot(0, 420, -5, 1, 100, 2);
        auto* legend = plot->getChart()->legend();
        legend->setVisible(true);
        legend->setAlignment(Qt::AlignTop);
        legend->setMarkerShape(QLegend::MarkerShapeFromSeries);
        auto* axisX = plot->getAxisX();
        auto* axisY = plot->getAxisY();
        axisX->setTitleText("Step");
        axisX->setLabelFormat("%d");
        axisY->setTitleText("Force/a.u.");
        axisY->setLabelFormat("10<sup>%d</sup>");
        {
            Vector<ScalarType> bfgs;
            PWscfOut pwout("./bfgs.out", 12);
            bfgs = pwout.makeTotalForces();
            const VectorType logNorm = ln(bfgs) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("BFGS(QE)");
        }
        {
            Vector<ScalarType> damp;
            PWscfOut pwout("./damp.out", 12);
            damp = pwout.makeTotalForces();
            const VectorType logNorm = ln(damp) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("Damp(QE)");
        }
        {
            Vector<ScalarType> fire;
            PWscfOut pwout("./fire.out", 12);
            fire = pwout.makeTotalForces();
            const VectorType logNorm = ln(fire) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("FIRE(QE)");
        }
        {
            Vector<ScalarType> fire;
            H5File h5f("FireQE.h5", H5File::ReadOnly);
            fire.read(h5f, "force");
            const VectorType logNorm = ln(fire) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("FIRE1(Physica)");
        }
        {
            Vector<ScalarType> fire;
            H5File h5f("FireQE1.h5", H5File::ReadOnly);
            fire.read(h5f, "force");
            const VectorType logNorm = ln(fire) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("FIRE2(Physica)");
        }
        plot->show();
    }
    return QApplication::exec();
}
