/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtCharts>
#include <Physica/Core/Math/Calculus/Differential.h>
#include "Physica/Gui/Plot.h"

using namespace Physica::Core;

namespace Physica::Gui {
    Plot::Plot(MultiScalar (*func)(const MultiScalar&), const MultiScalar& begin, const MultiScalar& end, QWidget* parent)
            : QtCharts::QChartView(parent), series(new QtCharts::QSplineSeries()) {
        setAttribute(Qt::WA_DeleteOnClose);

        MultiScalar maxStepSize = (end - begin) / BasicConst::getInstance().getPlotPoints();
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