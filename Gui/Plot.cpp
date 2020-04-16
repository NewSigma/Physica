/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Plot.h"
#include "Numerical.h"
#include <QtCharts>

Plot::Plot(Numerical (*func)(const Numerical&), const Numerical& begin, const Numerical& end, QWidget* parent)
    : QtCharts::QChartView(parent), series(new QtCharts::QSplineSeries()) {
    Numerical plotStepSize = (end - begin) / basicConst->getPlotPoints();
    Numerical point = Numerical(begin);
    for(unsigned long l = basicConst->getPlotPoints()[0]; l > 1; --l) {
        Numerical y = func(point);
        *series << QPointF(double(point), double(y));
        point += plotStepSize;
    }
    //Handle l = 1.
    Numerical y = func(point);
    *series << QPointF(double(point), double(y));

    auto chart = new QChart();
    chart->addSeries(series);
    chart->createDefaultAxes();

    setChart(chart);
    setRenderHint(QPainter::Antialiasing);
}