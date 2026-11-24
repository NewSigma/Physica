/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include <QtCharts/QChartView>
#include <QtCharts/QPolarChart>
#include <QtCharts/QSplineSeries>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/VectorImpl/RValueVector.h"

namespace Physica {
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
