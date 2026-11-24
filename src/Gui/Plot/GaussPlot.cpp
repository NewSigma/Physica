/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2023-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "Physica/Gui/Plot/GaussPlot.h"

using namespace Physica;

GaussPlot::GaussPlot(double maxX, double maxY, double deltaX, double deltaY, unsigned int numSigma)
        : Base(-maxX, maxX, 0, maxY, deltaX, deltaY) {
    using namespace Physica;
    using ScalarType = float64;
    using VectorType = VectorND<ScalarType>;

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
