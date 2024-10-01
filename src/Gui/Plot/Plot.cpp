/*
 * Copyright 2019-2024 Weibo He.
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
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;

namespace Physica::Gui {
    Plot::Plot(double minX, double maxX, double minY, double maxY, double deltaX, double deltaY, QWidget* parent)
            : ChartView(parent)
            , axisX(new QValueAxis())
            , axisY(new QValueAxis())
            , axisTop(new QValueAxis())
            , axisRight(new QValueAxis()) {
        setAttribute(Qt::WA_DeleteOnClose);

        setChart(new QChart());
        setRenderHint(QPainter::Antialiasing);
        setBackgroundRole(QPalette::Light);

        auto& chart = *Base::getChart();
        auto* legend = chart.legend();
        QFont font = legend->font();
        font.setPointSize(15);
        legend->setFont(font);
        legend->setAlignment(Qt::AlignTop);
        legend->setMarkerShape(QLegend::MarkerShapeFromSeries);
        chart.setTitleFont(font);
        chart.setBackgroundVisible(false);
        chart.setMargins(QMargins{});
        {
            axisX->setTickAnchor(0);
            axisX->setTickInterval(deltaX);
            axisX->setTickType(QValueAxis::TicksDynamic);
            axisX->setGridLineVisible(false);
            axisX->setMinorGridLineVisible(false);
            axisX->setLabelsFont(font);
            axisX->setRange(minX, maxX);
            axisX->setTitleFont(font);
            axisX->setLinePenColor(Qt::black);

            axisY->setTickAnchor(0);
            axisY->setTickInterval(deltaY);
            axisY->setTickType(QValueAxis::TicksDynamic);
            axisY->setGridLineVisible(false);
            axisY->setMinorGridLineVisible(false);
            axisY->setMinorTickCount(4);
            axisY->setLabelsFont(font);
            axisY->setRange(minY, maxY);
            axisY->setTitleFont(font);
            axisY->setLinePenColor(Qt::black);

            axisTop->setTickAnchor(0);
            axisTop->setTickInterval(deltaX);
            axisTop->setTickType(QValueAxis::TicksDynamic);
            axisTop->setLabelsVisible(false);
            axisTop->setGridLineVisible(false);
            axisTop->setMinorGridLineVisible(false);
            axisTop->setLabelsFont(font);
            axisTop->setRange(minX, maxX);
            axisTop->setTitleFont(font);
            axisTop->setLinePenColor(Qt::black);

            axisRight->setTickAnchor(0);
            axisRight->setTickInterval(deltaY);
            axisRight->setTickType(QValueAxis::TicksDynamic);
            axisRight->setLabelsVisible(false);
            axisRight->setGridLineVisible(false);
            axisRight->setMinorGridLineVisible(false);
            axisRight->setMinorTickCount(4);
            axisRight->setLabelsFont(font);
            axisRight->setRange(minY, maxY);
            axisRight->setTitleFont(font);
            axisRight->setLinePenColor(Qt::black);

            chart.addAxis(axisX, Qt::AlignBottom);
            chart.addAxis(axisY, Qt::AlignLeft);
            chart.addAxis(axisTop, Qt::AlignTop);
            chart.addAxis(axisRight, Qt::AlignRight);
        }
    }

    QScatterSeries& Plot::label(double x, double y, QString text) {
        using VectorType = Vector<float64, 1>;
        auto& result = scatter(VectorType{x}, VectorType{y});
        result.setPointLabelsVisible(true);
        result.setPointLabelsFormat(QPointLabelFormat(std::move(text)));
        result.setMarkerSize(0);
        return result;
    }
}