/*
 * Copyright 2022-2026 Weibo He.
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
module;

#include <QtCharts/QChartView>
#include <QtCharts/QPolarChart>
#include <QtCharts/QSplineSeries>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"

export module Physica.Gui.PolarPlot;

export namespace Physica {
    class PHYSICA_API PolarPlot : public QChartView {
        using Base = QChartView;
    public:
        PolarPlot(QWidget* parent = nullptr);

        QSplineSeries& spline(const Vector auto& theta, const Vector auto& rho);
        /* Getters */
        [[nodiscard]] QPolarChart* chart() { return static_cast<QPolarChart*>(Base::chart()); }
    };

    QSplineSeries& PolarPlot::spline(const Vector auto& theta, const Vector auto& rho) {
        assert(rho.getLength() == theta.getLength());
        QSplineSeries* series = new QSplineSeries();
        for (size_t i = 0; i < rho.getLength(); ++i)
            *series << QPointF(double(theta.calc(i)), double(rho.calc(i)));
        chart()->addSeries(series);

        update();
        return *series;
    }
}
