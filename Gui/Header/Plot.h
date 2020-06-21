/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PLOT_H
#define PHYSICA_PLOT_H

#include <QChartView>
#include <QtCharts/QSplineSeries>

namespace Physica::Core {
    class Scalar;
}
using Physica::Core::Scalar;

class Plot : public QtCharts::QChartView {
    QtCharts::QSplineSeries* series;
public:
    Plot(Scalar (*func)(const Scalar&), const Scalar& begin, const Scalar& end, QWidget* parent = nullptr);
};

#endif
