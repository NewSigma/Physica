/*
 * Copyright 2019-2021 WeiBo He.
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
#pragma once

#include <QtCharts>
#include <QChartView>
#include <QtCharts/QSplineSeries>
#include <Physica/Core/MultiPrecision/Scalar.h>
#include <Physica/Utils/Container/Array/Array.h>

using Physica::Core::MultiScalar;

namespace Physica::Gui {
    class Plot : public QtCharts::QChartView {
        QtCharts::QSplineSeries* series;
    public:
        Plot(MultiScalar (*func)(const MultiScalar&), const MultiScalar& begin
                , const MultiScalar& end, QWidget* parent = nullptr);
        template<class T, size_t Length, size_t Capacity>
        Plot(const Utils::Array<T, Length, Capacity>& x, const Utils::Array<T, Length, Capacity>& y, QWidget* parent = nullptr);
    };

    template<class T, size_t Length, size_t Capacity>
    Plot::Plot(const Utils::Array<T, Length, Capacity>& x, const Utils::Array<T, Length, Capacity>& y, QWidget* parent)
            : QtCharts::QChartView(parent), series(new QtCharts::QSplineSeries()) {
        setAttribute(Qt::WA_DeleteOnClose);
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x[i]), double(y[i]));

        auto chart = new QChart();
        chart->addSeries(series);
        chart->createDefaultAxes();

        setChart(chart);
        setRenderHint(QPainter::Antialiasing);
    }
}
