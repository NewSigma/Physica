/*
 * Copyright 2023-2026 Weibo He.
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
#include "Physica/Core/Physics/SolidState/QE/PWscfOut.h"

import Physica.Gui.Plot;

using namespace Physica;
using ScalarType = float64;
using VectorType = VectorND<ScalarType>;

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
            VectorND<ScalarType> bfgs;
            PWscfOut pwout("./bfgs.out", 12);
            bfgs = pwout.makeTotalForces();
            const VectorType logNorm = ln(bfgs) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("BFGS(QE)");
        }
        {
            VectorND<ScalarType> damp;
            PWscfOut pwout("./damp.out", 12);
            damp = pwout.makeTotalForces();
            const VectorType logNorm = ln(damp) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("Damp(QE)");
        }
        {
            VectorND<ScalarType> fire;
            PWscfOut pwout("./fire.out", 12);
            fire = pwout.makeTotalForces();
            const VectorType logNorm = ln(fire) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("FIRE(QE)");
        }
        {
            VectorND<ScalarType> fire;
            auto h5f = H5File::open("FireQE.h5", H5File::ReadOnly);
            fire.read(h5f, "force");
            const VectorType logNorm = ln(fire) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("FIRE1(Physica)");
        }
        {
            VectorND<ScalarType> fire;
            auto h5f = H5File::open("FireQE1.h5", H5File::ReadOnly);
            fire.read(h5f, "force");
            const VectorType logNorm = ln(fire) * reciprocal(ln(ScalarType(10)));
            plot->line(logNorm).setName("FIRE2(Physica)");
        }
        plot->show();
    }
    return QApplication::exec();
}
