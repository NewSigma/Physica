/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2020-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <iostream>
#include <QtWidgets/QApplication>
#include "Physica/Core/Math/Geometry/Point.h"
#include "Physica/Gui/Plot/DensityPlot.h"

int main(int argc, char** argv) {
    using namespace Physica;

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
