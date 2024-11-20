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
#include <QSvgGenerator>
#include "Physica/Gui/Plot/Plot.h"

using namespace Physica::Core;

namespace Physica::Gui {
    Plot::Plot(QWidget* parent) : ChartView(parent), axisX(new QValueAxis()), axisY(new QValueAxis()), axisTop(new QValueAxis()), axisRight(new QValueAxis()) {
        setAttribute(Qt::WA_DeleteOnClose);

        setChart(new QChart());
        setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing | QPainter::SmoothPixmapTransform);
        setBackgroundRole(QPalette::Light);

        auto& legend = getLegend();
        legend.setAlignment(Qt::AlignTop);
        legend.setMarkerShape(QLegend::MarkerShapeFromSeries);
        legend.hide();

        auto& chart = *Base::getChart();
        chart.setBackgroundVisible(false);
        chart.setMargins(QMargins{});
        {
            axisX->setTickAnchor(0);
            axisX->setTickType(QValueAxis::TicksDynamic);
            axisX->setGridLineVisible(false);
            axisX->setMinorGridLineVisible(false);
            axisX->setLinePenColor(Qt::black);

            axisY->setTickAnchor(0);
            axisY->setTickType(QValueAxis::TicksDynamic);
            axisY->setGridLineVisible(false);
            axisY->setMinorGridLineVisible(false);
            axisY->setLinePenColor(Qt::black);

            axisTop->setTickAnchor(0);
            axisTop->setTickType(QValueAxis::TicksDynamic);
            axisTop->setLabelsVisible(false);
            axisTop->setGridLineVisible(false);
            axisTop->setMinorGridLineVisible(false);
            axisTop->setLinePenColor(Qt::black);

            axisRight->setTickAnchor(0);
            axisRight->setTickType(QValueAxis::TicksDynamic);
            axisRight->setLabelsVisible(false);
            axisRight->setGridLineVisible(false);
            axisRight->setMinorGridLineVisible(false);
            axisRight->setLinePenColor(Qt::black);

            setTickDirection(QAbstractAxis::Inner);
            chart.addAxis(axisX, Qt::AlignBottom);
            chart.addAxis(axisY, Qt::AlignLeft);
            chart.addAxis(axisTop, Qt::AlignTop);
            chart.addAxis(axisRight, Qt::AlignRight);
        }

        QFont font{};
        font.setPointSize(15);
        setFont(font);
    }

    Plot::Plot(double minX, double maxX, double minY, double maxY, double deltaX, double deltaY, QWidget* parent) : Plot(parent) {
        setBox(minX, maxX, minY, maxY, deltaX, deltaY);
    }

    QScatterSeries& Plot::label(double x, double y, QString text) {
        using Vec = Vector1D<float64>;
        auto& result = scatter(Vec{x}, Vec{y});
        result.setPointLabelsVisible(true);
        result.setPointLabelsFormat(QPointLabelFormat(std::move(text)));
        result.setMarkerSize(0);
        return result;
    }

    QImage Plot::toImage(int width, int height) {
        auto image = QImage(QSize(width, height), QImage::Format_ARGB32);

        QPainter painter{};
        painter.begin(&image);
        painter.setRenderHints(renderHints());
        painter.fillRect(image.rect(), Qt::white);

        auto& chart = *getChart();
        const auto size1 = image.size() * 0.99;
        chart.setMinimumSize(size1);
        chart.setMaximumSize(size1);

        render(&painter, image.rect(), chart.rect().toRect());
        painter.end();

        chart.setMaximumSize(0, 0);
        chart.setMaximumSize(QWIDGETSIZE_MAX, QWIDGETSIZE_MAX);
        return image;
    }

    void Plot::toSvg(const char* path, int width, int height, int resolution) {
        QSvgGenerator gen{};
        gen.setFileName(path);
        gen.setTitle("SVG Plot");
        gen.setDescription("Created by Physica(https://gitee.com/newsigma/Physica)");
        gen.setSize(QSize(width, height));
        gen.setViewBox(QRectF(QPointF{}, gen.size()));
        if (resolution > 0)
            gen.setResolution(resolution);

        QPainter painter{};
        painter.begin(&gen);
        painter.setRenderHints(renderHints());
        painter.fillRect(gen.viewBox(), Qt::white);

        auto& chart = *getChart();
        const auto size1 = gen.size() * 0.99;
        chart.setMinimumSize(size1);
        chart.setMaximumSize(size1);

        render(&painter);
        painter.end();

        chart.setMaximumSize(0, 0);
        chart.setMaximumSize(QWIDGETSIZE_MAX, QWIDGETSIZE_MAX);
    }

    void Plot::setBox(double minX, double maxX, double minY, double maxY, double deltaX, double deltaY) {
        axisX->setTickInterval(deltaX);
        axisX->setRange(minX, maxX);
        axisY->setTickInterval(deltaY);
        axisY->setRange(minY, maxY);
        axisTop->setTickInterval(deltaX);
        axisTop->setRange(minX, maxX);
        axisRight->setTickInterval(deltaY);
        axisRight->setRange(minY, maxY);
    }

    void Plot::setTickDirection(QAbstractAxis::TickDirection d) {
        axisX->setTickDirection(d);
        axisY->setTickDirection(d);
        axisTop->setTickDirection(d);
        axisRight->setTickDirection(d);
    }

    void Plot::setFont(QFont font) {
        Base::setFont(font);
        getLegend().setFont(font);
        getChart()->setTitleFont(font);

        axisX->setLabelsFont(font);
        axisX->setTitleFont(font);
        axisY->setLabelsFont(font);
        axisY->setTitleFont(font);
        axisTop->setLabelsFont(font);
        axisTop->setTitleFont(font);
        axisRight->setLabelsFont(font);
        axisRight->setTitleFont(font);
    }

    void Plot::setFontSize(int size) {
        auto f = Base::font();
        f.setPointSize(size);
        setFont(f);
    }
}
