/*
 * Copyright 2023 Weibo He.
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
#include "Physica/Gui/Plot/GaussPlot.h"

namespace Physica::Gui {
    GaussPlot::GaussPlot(double maxX, double maxY, double deltaX, double deltaY, unsigned int numSigma)
            : Base(-maxX, maxX, 0, maxY, deltaX, deltaY) {
        using namespace Physica::Core;
        using ScalarType = float64;
        using VectorType = Vector<ScalarType>;

        auto& line = Base::line(VectorType{-maxX, maxX}, VectorType{0, 0});
        auto pen = line.pen();
        pen.setColor(Qt::black);
        pen.setStyle(Qt::DashLine);
        line.setPen(pen);

        auto& line1 = Base::line(VectorType{0, 0}, VectorType{0, maxY});
        line1.setPen(pen);

        pen.setColor(Qt::blue);
        auto& line2 = Base::line(VectorType{0, maxX}, VectorType{0, maxX / numSigma});
        line2.setPen(pen);

        auto& line3 = Base::line(VectorType{-maxX, 0}, VectorType{maxX / numSigma, 0});
        line3.setPen(pen);
    }
}
