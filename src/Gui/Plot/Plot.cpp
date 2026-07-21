/*
 * Copyright 2019-2026 Weibo He.
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

#include <QWidget>
#include <QPainter>
#include <QSvgGenerator>
#include <QtCharts/QChart>
#include <QtCharts/QValueAxis>
#include <QtCharts/QScatterSeries>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"

module Physica.Gui.Plot;

using namespace Physica;

Plot::Plot(QWidget* parent)
        : ChartView(parent), axisX(new QValueAxis()), axisY(new QValueAxis()), axisTop(new QValueAxis()), axisRight(new QValueAxis()) {
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

Plot::Plot(float64 minX, float64 maxX, float64 minY, float64 maxY, float64 deltaX, float64 deltaY, QWidget* parent)
        : Plot(parent) {
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

void Plot::setBox(float64 minX, float64 maxX, float64 minY, float64 maxY, float64 deltaX, float64 deltaY) {
    assert(deltaX.isPositive());
    assert(deltaY.isPositive());
    axisX->setTickInterval((double)deltaX);
    axisX->setRange((double)minX, (double)maxX);
    axisY->setTickInterval((double)deltaY);
    axisY->setRange((double)minY, (double)maxY);
    axisTop->setTickInterval((double)deltaX);
    axisTop->setRange((double)minX, (double)maxX);
    axisRight->setTickInterval((double)deltaY);
    axisRight->setRange((double)minY, (double)maxY);
}

void Plot::setAxisX(QValueAxis* axis) {
    auto* chart = Base::getChart();
    chart->removeAxis(axisX);
    chart->addAxis(axis, Qt::AlignBottom);
    axisX = axis;
}

void Plot::setAxisY(QValueAxis* axis) {
    auto* chart = Base::getChart();
    chart->removeAxis(axisY);
    chart->addAxis(axis, Qt::AlignLeft);
    axisY = axis;
}

void Plot::setAxisTop(QValueAxis* axis) {
    auto* chart = Base::getChart();
    chart->removeAxis(axisTop);
    chart->addAxis(axis, Qt::AlignTop);
    axisTop = axis;
}

void Plot::setAxisRight(QValueAxis* axis) {
    auto* chart = Base::getChart();
    chart->removeAxis(axisRight);
    chart->addAxis(axis, Qt::AlignRight);
    axisRight = axis;
}

void Plot::setMinX(double value) noexcept {
    axisX->setMin(value);
    axisTop->setMin(value);
}

void Plot::setMaxX(double value) noexcept {
    axisX->setMax(value);
    axisTop->setMax(value);
}

void Plot::setRangeX(double minX, double maxX) {
    setMinX(minX);
    setMaxX(maxX);
}

void Plot::setMinY(double value) noexcept {
    axisY->setMin(value);
    axisRight->setMin(value);
}

void Plot::setMaxY(double value) noexcept {
    axisY->setMax(value);
    axisRight->setMax(value);
}

void Plot::setRangeY(double minY, double maxY) {
    setMinY(minY);
    setMaxY(maxY);
}

void Plot::setDeltaX(double value) noexcept {
    axisX->setTickInterval(value);
    axisTop->setTickInterval(value);
}

void Plot::setDeltaY(double value) noexcept {
    axisY->setTickInterval(value);
    axisRight->setTickInterval(value);
}

void Plot::setTickDirection(QAbstractAxis::TickDirection d) {
    axisX->setTickDirection(d);
    axisY->setTickDirection(d);
    axisTop->setTickDirection(d);
    axisRight->setTickDirection(d);
}

void Plot::setFont(const QFont& font) {
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
