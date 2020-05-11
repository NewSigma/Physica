/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PLOT_H
#define PHYSICA_PLOT_H

#include <QChartView>
#include <QtCharts/QSplineSeries>

namespace Physica::Core {
    class Numerical;
}
using Physica::Core::Numerical;

class Plot : public QtCharts::QChartView {
    QtCharts::QSplineSeries* series;
public:
    Plot(Numerical (*func)(const Numerical&), const Numerical& begin, const Numerical& end, QWidget* parent = nullptr);
};

#endif
