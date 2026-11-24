/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <QApplication>
#include "Physica/Core/Physics/SolidState/QE/PWscfOut.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;

int main(int argc, char** argv) {
    QApplication app(argc, argv);
    {
        Plot* plot = new Plot(0, 200, -7, 1, 50, 2);
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
            VectorND<ScalarType> bfgs;
            PWscfOut pwout("./bfgs1.out", 12);
            bfgs = pwout.makeTotalForces();
            const VectorType logNorm = ln(bfgs) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("BFGS(QE)");
        }
        {
            VectorND<ScalarType> damp;
            PWscfOut pwout("./damp1.out", 12);
            damp = pwout.makeTotalForces();
            const VectorType logNorm = ln(damp) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("Damp(QE)");
        }
        {
            VectorND<ScalarType> fire;
            auto h5f = H5File::open("FireQE2.h5", H5File::ReadOnly);
            fire.read(h5f, "force");
            const VectorType logNorm = ln(fire) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("FIRE(Physica)");
        }
        plot->show();
    }
    return QApplication::exec();
}
