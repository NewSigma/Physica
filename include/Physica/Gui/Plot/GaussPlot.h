/*
 * Copyright 2023 WeiBo He.
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

#include <QtCharts/QValueAxis>
#include "Plot.h"

namespace Physica::Gui {
    class GaussPlot : public Plot {
        using Base = Plot;

        double minX;
        double maxX;
        double minY;
        double maxY;
        double deltaX;
        double deltaY;
        QValueAxis* axisX;
        QValueAxis* axisY;
        QValueAxis* axisTop;
        QValueAxis* axisRight;
    public:
        GaussPlot(double minX_, double maxX_, double minY_, double maxY_, double deltaX_, double deltaY_, unsigned int numSigma);
        /* Getters */
        [[nodiscard]] QValueAxis* getAxisX() const noexcept { return axisX; }
        [[nodiscard]] QValueAxis* getAxisY() const noexcept { return axisY; }
        [[nodiscard]] QValueAxis* getAxisTop() const noexcept { return axisTop; }
        [[nodiscard]] QValueAxis* getAxisRight() const noexcept { return axisRight; }
    };
}
