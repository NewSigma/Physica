/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

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
#include <QtCharts>
#include <Physica/Core/Math/Calculus/Differential.h>
#include "Physica/Gui/Plot.h"

using namespace Physica::Core;

namespace Physica::Gui {
    Plot::Plot(MultiScalar (*func)(const MultiScalar&), const MultiScalar& begin, const MultiScalar& end, QWidget* parent)
            : QtCharts::QChartView(parent), series(new QtCharts::QSplineSeries()) {
        setAttribute(Qt::WA_DeleteOnClose);

        MultiScalar maxStepSize = (end - begin) / BasicConst::getInstance().plotPoints;
        MultiScalar point = MultiScalar(begin);
        do {
            MultiScalar y = func(point);
            *series << QPointF(double(point), double(y));
            /*
            //Changeable step size depending on current derivative.
            MultiScalar derivative = D_right(type, point);
            std::cout << derivative << point << std::endl;
            if(derivative.getPower() > basicConst->MaxPower) {
                qWarning("StepSize is too small, suspect encountered a singularity. Ignoring...");
                derivative = basicConst->get_0();
            }
            point += maxStepSize / (basicConst->get_1() + derivative * derivative);
             */
            point += maxStepSize;
            point.clearA();
        } while(point < end);
        //Handle l = 1.
        MultiScalar y = func(point);
        *series << QPointF(double(point), double(y));

        auto chart = new QChart();
        chart->addSeries(series);
        chart->createDefaultAxes();

        setChart(chart);
        setRenderHint(QPainter::Antialiasing);
    }
}