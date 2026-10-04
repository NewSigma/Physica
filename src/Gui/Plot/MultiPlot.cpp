/*
 * Copyright 2024-2026 Weibo He.
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
module;

#include <QSvgGenerator>
#include <QGridLayout>

module Physica.Gui.MultiPlot;

using namespace Physica;

MultiPlot::MultiPlot(int row, int col, QWidget* parent)
        : QWidget(parent), layout(new QGridLayout(this)) {
    QPalette p{};
    p.setColor(QPalette::Window, Qt::white);
    setPalette(p);
    setAutoFillBackground(true);

    for (int r = 0; r < row; ++r)
        for (int c = 0; c < col; ++c)
            layout->addWidget(new Plot(), r, c);
}

MultiPlot::MultiPlot(int row, int col, double minX, double maxX, double minY, double maxY, double deltaX, double deltaY, QWidget* parent)
        : MultiPlot(row, col, parent) {
    setBox(minX, maxX, minY, maxY, deltaX, deltaY);
}

Plot& MultiPlot::operator[](int row, int col) {
    assert(0 <= row && row < getRow());
    assert(0 <= col && col < getCol());
    return static_cast<Plot&>(*layout->itemAtPosition(row, col)->widget());
}

void MultiPlot::toSvg(const char* path, int resolution) {
    QSvgGenerator gen{};
    gen.setFileName(path);
    gen.setTitle("SVG Plot");
    gen.setDescription("Created by Physica(https://gitee.com/newsigma/Physica)");
    gen.setSize(Base::size());
    gen.setViewBox(QRectF(QPointF{}, gen.size()));
    if (resolution > 0)
        gen.setResolution(resolution);

    QPainter painter{};
    painter.begin(&gen);
    painter.setRenderHints(operator[](0, 0).renderHints());
    painter.fillRect(gen.viewBox(), Qt::white);
    render(&painter);
    painter.end();
}
