/*
 * Copyright 2025-2026 Weibo He.
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
#include <QApplication>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

import Physica.Gui.Plot;

using namespace Physica;

int main(int argc, char** argv) {
    QApplication app(argc, argv); // Initialize Qt

    Plot* plot = new Plot(-5, 5, -1.1, 1.1, 2, 0.5);
    plot->getAxisX()->setLabelFormat("%d");
    plot->getAxisY()->setLabelFormat("%.1f");
    plot->getAxisX()->setTitleText("x");
    plot->getAxisY()->setTitleText("y");

    auto x = VectorND<float64>::linspace(-5, 5, 100);
    plot->spline(x, tanh(x));

    auto& line = plot->line(Vector2D<float64>{-5, 5}, Vector2D<float64>{0, 0}); // Line from (-5, 0) to (5, 0)
    auto p = line.pen();
    p.setColor(Qt::black);
    p.setStyle(Qt::DashLine);
    line.setPen(p);

    auto& l = plot->label(3, 0.7, "y = tanh(x)");
    auto font = l.pointLabelsFont();
    font.setPointSize(20);
    l.setPointLabelsFont(font);

    plot->show();
    return QApplication::exec(); // Wait until user close plot
}
