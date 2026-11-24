/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include <QSvgGenerator>
#include "Physica/Gui/Plot/MultiPlot.h"

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
