/*
 * Copyright 2024 WeiBo He.
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
#include "Physica/Gui/Plot/PlotImpl/ChartView.h"

namespace Physica::Gui {
    ChartView::ChartView(QChart* chart_, QWidget* parent)
            : QGraphicsView(parent)
            , chart(chart_)
            , rubberBand(*new QRubberBand(QRubberBand::Rectangle, this)) {
        rubberBand.setEnabled(true);

        setFrameShape(QFrame::NoFrame);
        setBackgroundRole(QPalette::Window);
        setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        setScene(new QGraphicsScene(this));
        setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        if (!chart)
            chart = new QChart();
        scene()->addItem(chart);
    }

    void ChartView::resizeEvent(QResizeEvent* event) {
        QGraphicsView::resizeEvent(event);
        resize();
    }

    void ChartView::setChart(QChart* chart_) {
        assert(chart_ != nullptr);
        if (chart == chart_)
            return;
        if (chart)
            scene()->removeItem(chart);

        chart = chart_;
        scene()->addItem(chart);
        resize();
    }

    void ChartView::mousePressEvent(QMouseEvent* event) {
        const QRectF plotArea = chart->plotArea();
        if (rubberBand.isEnabled() && event->button() == Qt::LeftButton && plotArea.contains(event->pos())) {
            rubberBandOrigin = event->pos();
            rubberBand.setGeometry(QRect(rubberBandOrigin, QSize()));
            rubberBand.show();
            event->accept();
        }
        else
            Base::mousePressEvent(event);
    }

    void ChartView::mouseMoveEvent(QMouseEvent* event) {
        if (rubberBand.isVisible()) {
            const int x = rubberBandOrigin.x();
            const int y = rubberBandOrigin.y();
            const int width = event->pos().x() - x;
            const int height = event->pos().y() - y;
            rubberBand.setGeometry(QRect(x, y, width, height).normalized());
        }
        else
            QGraphicsView::mouseMoveEvent(event);
    }

    void ChartView::mouseReleaseEvent(QMouseEvent* event) {
        if (rubberBand.isVisible()) {
            if (event->button() == Qt::LeftButton) {
                rubberBand.hide();
                QRectF rect = rubberBand.geometry();
                chart->zoomIn(rect);
                event->accept();
            }
        }
        else if (event->button() == Qt::RightButton) {
            chart->zoomReset();
            event->accept();
        }
        else
            QGraphicsView::mouseReleaseEvent(event);
    }

    void ChartView::resize() {
        QSize chartSize = size();
        const double sinA = qAbs(transform().m21());
        const bool isRotated = sinA != 1.0;
        if (isRotated) {
            const bool isNotDegree90 = sinA != 0.0;
            if (isNotDegree90) {
                const double cosA = qAbs(transform().m11());
                double minDimension = qMin(size().width(), size().height());
                double h = (minDimension - (minDimension / ((sinA / cosA) + 1.0))) / sinA;
                chartSize.setHeight(h);
                chartSize.setWidth(h);
            }
        }
        else {
            chartSize.setHeight(size().width());
            chartSize.setWidth(size().height());
        }

        chart->resize(chartSize);
        setMinimumSize(chart->minimumSize().toSize().expandedTo(minimumSize()));
        setMaximumSize(maximumSize().boundedTo(chart->maximumSize().toSize()));
        setSceneRect(chart->geometry());
    }
}
