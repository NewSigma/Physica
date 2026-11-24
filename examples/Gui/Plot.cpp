/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <QApplication>
#include "Physica/Gui/Plot/Plot.h"

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
