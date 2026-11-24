/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2022-2025 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Gui/Plot/PolarPlot.h"

using namespace Physica;

PolarPlot::PolarPlot(QWidget* parent) : QChartView(parent) {
    setAttribute(Qt::WA_DeleteOnClose);

    setChart(new QPolarChart());
    setRenderHint(QPainter::Antialiasing);
}
