/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PLOT_H
#define PHYSICA_PLOT_H

#include <QChartView>
#include <QtCharts/QSplineSeries>
#include <Physica/Core/MultiPrecition/Scalar.h>

using Physica::Core::MultiScalar;

namespace Physica::Gui {
    class Plot : public QtCharts::QChartView {
        QtCharts::QSplineSeries* series;
    public:
        Plot(MultiScalar (*func)(const MultiScalar&), const MultiScalar& begin
                , const MultiScalar& end, QWidget* parent = nullptr);
    };
}

#endif
