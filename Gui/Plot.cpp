/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Plot.h"
#include "Core/Header/Scalar.h"
#include <QtCharts>
#include <Core/Header/Differential.h>

using namespace Physica::Core;

Plot::Plot(Scalar (*func)(const Scalar&), const Scalar& begin, const Scalar& end, QWidget* parent)
    : QtCharts::QChartView(parent), series(new QtCharts::QSplineSeries()) {
    setAttribute(Qt::WA_DeleteOnClose);

    Scalar maxStepSize = (end - begin) / BasicConst::getInstance().getPlotPoints();
    Scalar point = Scalar(begin);
    do {
        Scalar y = func(point);
        *series << QPointF(double(point), double(y));
        /*
        //Changeable step size depending on current derivative.
        Scalar derivative = D_right(func, point);
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
    Scalar y = func(point);
    *series << QPointF(double(point), double(y));

    auto chart = new QChart();
    chart->addSeries(series);
    chart->createDefaultAxes();

    setChart(chart);
    setRenderHint(QPainter::Antialiasing);
}