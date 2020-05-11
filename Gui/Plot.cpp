/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Plot.h"
#include "Core/Header/Numerical.h"
#include <QtCharts>
#include <Core/Header/Differential.h>

using namespace Physica::Core;

Plot::Plot(Numerical (*func)(const Numerical&), const Numerical& begin, const Numerical& end, QWidget* parent)
    : QtCharts::QChartView(parent), series(new QtCharts::QSplineSeries()) {
    setAttribute(Qt::WA_DeleteOnClose);

    Numerical maxStepSize = (end - begin) / BasicConst::getInstance().getPlotPoints();
    Numerical point = Numerical(begin);
    do {
        Numerical y = func(point);
        *series << QPointF(double(point), double(y));
        /*
        //Changeable step size depending on current derivative.
        Numerical derivative = D_right(func, point);
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
    Numerical y = func(point);
    *series << QPointF(double(point), double(y));

    auto chart = new QChart();
    chart->addSeries(series);
    chart->createDefaultAxes();

    setChart(chart);
    setRenderHint(QPainter::Antialiasing);
}