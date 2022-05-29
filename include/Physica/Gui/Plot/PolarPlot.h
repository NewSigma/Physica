/*
 * Copyright 2022 WeiBo He.
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

#include <QtCharts/QChartView>
#include <QtCharts/QPolarChart>
#include <QtCharts/QSplineSeries>
#include <Physica/Core/Math/Algebra/LinearAlgebra/Vector/RValueVector.h>

namespace Physica::Gui {
    class PolarPlot : public QChartView {
        using Base = QChartView;
    public:
        PolarPlot(QWidget* parent = nullptr);

        template<class VectorType1, class VectorType2>
        QSplineSeries& spline(const Core::RValueVector<VectorType1>& theta, const Core::RValueVector<VectorType2>& rho);
        /* Getters */
        [[nodiscard]] QPolarChart* chart() { return static_cast<QPolarChart*>(Base::chart()); }
    };

    template<class VectorType1, class VectorType2>
    QSplineSeries& PolarPlot::spline(const Core::RValueVector<VectorType1>& theta, const Core::RValueVector<VectorType2>& rho) {
        assert(rho.getLength() == theta.getLength());
        QSplineSeries* series = new QSplineSeries();
        for (size_t i = 0; i < rho.getLength(); ++i)
            *series << QPointF(double(theta.calc(i)), double(rho.calc(i)));
        chart()->addSeries(series);

        update();
        return *series;
    }
}
