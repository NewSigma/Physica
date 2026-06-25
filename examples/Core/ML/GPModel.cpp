/*
 * Copyright 2026 Weibo He.
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
#include "Physica/Core/ML/GPModel.h"
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica;
using T = float64;
using RandomSource = Random<>;
/**
 * Gaussian Process Regression for noisy sine function prediction
 */
int main(int argc, char** argv) {
    GPModel<T> gp(1);
    const auto x0 = VectorND<T>::linspace(-MathConst<T>::pi, MathConst<T>::pi, 8);
    const auto y0 = VectorND<T>(sin(x0) + VectorND<T>::random_normal<RandomSource>(x0.getLength()) * T(1E-3));
    gp.regression(x0.transpose(), y0, T(1E-3));

    QApplication app(argc, argv);
    Plot* plot = new Plot(-MathConst<T>::pi, MathConst<T>::pi, -1.4, 1.4, 1, 0.5);
    plot->getAxisX()->setLabelFormat("%d");
    plot->getAxisY()->setLabelFormat("%.1f");
    plot->getAxisX()->setTitleText("x");
    plot->getAxisY()->setTitleText("y");
    {
        const auto x = VectorND<T>::linspace(-MathConst<T>::pi, MathConst<T>::pi, 200);
        auto& l = plot->line(x, sin(x));
        QPen pen = l.pen();
        pen.setColor(Qt::black);
        pen.setStyle(Qt::DashLine);
        l.setPen(pen);
    }
    {
        auto& s = plot->scatter(x0, y0);
        s.setMarkerSize(8);
        s.setBorderColor(QColor(204, 51, 51));
        s.setColor(Qt::transparent);
        s.setMarkerShape(s.MarkerShapeTriangle);

        const auto x = VectorND<T>::linspace(-MathConst<T>::pi, MathConst<T>::pi, 100);
        const auto pairs = Array<Vector2D<T>>::generate(x.getLength(), [&](size_t i) { return gp.predict({x[i]}); });
        const auto y = VectorND<T>::generate(x.getLength(), [&](size_t i) { return pairs[i][0]; });
        const auto dy = VectorND<T>::generate(x.getLength(), [&](size_t i) { return pairs[i][1]; });
        auto& l = plot->line(x, y);
        l.setColor(QColor(51, 102, 204));
        
        auto& area = plot->area_center(x, y, dy);
        auto color = l.color();
        color.setAlpha(60);
        area.setColor(color);
        area.setBorderColor(color);
    }
    plot->show();
    return QApplication::exec();
}
