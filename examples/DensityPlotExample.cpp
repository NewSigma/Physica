/*
 * Copyright 2020 WeiBo He.
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
#include <iostream>
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Geometry/Point.h"
#include "Physica/Gui/Plot/DensityPlot.h"

int main(int argc, char** argv) {
    using namespace Physica::Core;
    using namespace Physica::Gui;

    const int size = 1000;
    const int step_size = 10;
    const int radius = 14;
    QApplication app(argc, argv);
    DensityPlot plot(size, size);
    for(size_t i = 0; i < size; i += step_size) {
        for(size_t j = 0; j < size; j += step_size)
            plot.appendPoint(i, j, radius, uchar((double(i) * i + double(j) * j) / 2000000.0 * 255.0));
    }
    plot.show();
    return QApplication::exec();
}