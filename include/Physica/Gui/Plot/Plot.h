/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

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
#ifndef PHYSICA_PLOT_H
#define PHYSICA_PLOT_H

#include <QChartView>
#include <QtCharts/QSplineSeries>
#include <Physica/Core/MultiPrecision/Scalar.h>

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
