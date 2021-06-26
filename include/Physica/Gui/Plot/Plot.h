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
#include <QtCharts/QScatterSeries>
#include <Physica/Core/MultiPrecision/Scalar.h>
#include <Physica/Utils/Container/Array/Array.h>

using Physica::Core::MultiScalar;

namespace Physica::Gui {
    class Plot : public QtCharts::QChartView {
    public:
        Plot(QWidget* parent = nullptr);
        Plot(MultiScalar (*func)(const MultiScalar&), const MultiScalar& begin
                , const MultiScalar& end, QWidget* parent = nullptr);
        /* Operations */
        template<class Array>
        void spline(const Array& x, const Array& y);
        template<class Array>
        void scatter(const Array& x, const Array& y);
    };

    template<class Array>
    void Plot::spline(const Array& x, const Array& y) {
        QtCharts::QSplineSeries* series = new QtCharts::QSplineSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x[i]), double(y[i]));
        chart()->addSeries(series);
        chart()->createDefaultAxes();

        update();
    }

    template<class Array>
    void Plot::scatter(const Array& x, const Array& y) {
        QtCharts::QScatterSeries* series = new QtCharts::QScatterSeries();
        for (size_t i = 0; i < x.getLength(); ++i)
            *series << QPointF(double(x[i]), double(y[i]));
        chart()->addSeries(series);
        chart()->createDefaultAxes();

        update();
    }
}
