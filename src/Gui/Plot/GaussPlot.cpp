/*
 * Copyright 2023 WeiBo He.
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
#include "Physica/Gui/Plot/GaussPlot.h"

namespace Physica::Gui {
    GaussPlot::GaussPlot(double minX_, double maxX_, double minY_, double maxY_, double deltaX_, double deltaY_, unsigned int numSigma)
            : minX(minX_)
            , maxX(maxX_)
            , minY(minY_)
            , maxY(maxY_)
            , deltaX(deltaX_)
            , deltaY(deltaY_)
            , axisX(new QValueAxis())
            , axisY(new QValueAxis())
            , axisTop(new QValueAxis())
            , axisRight(new QValueAxis()) {
        using namespace Physica::Core;
        using ScalarType = Scalar<Double, false>;
        using VectorType = Vector<ScalarType>;
        auto& chart = *Base::chart();
        chart.legend()->setVisible(false);
        {
            QFont font = axisX->labelsFont();
            font.setPointSize(15);
            axisX->setTickAnchor(0);
            axisX->setTickInterval(deltaX);
            axisX->setTickType(QValueAxis::TicksDynamic);
            axisX->setMinorGridLineVisible(false);
            axisX->setLinePenColor(Qt::black);
            axisX->setGridLineVisible(false);
            axisX->setLabelsFont(font);
            axisX->setRange(minX, maxX);
            axisX->setTitleFont(font);

            axisY->setTickAnchor(0);
            axisY->setTickInterval(deltaY);
            axisY->setTickType(QValueAxis::TicksDynamic);
            axisY->setMinorGridLineVisible(false);
            axisY->setMinorTickCount(4);
            axisY->setLinePenColor(Qt::black);
            axisY->setGridLineVisible(false);
            axisY->setMinorGridLineVisible(false);
            axisY->setLabelsFont(font);
            axisY->setRange(minY, maxY);
            axisY->setTitleFont(font);

            axisTop->setTickAnchor(0);
            axisTop->setTickInterval(deltaX);
            axisTop->setTickType(QValueAxis::TicksDynamic);
            axisTop->setLabelsVisible(false);
            axisTop->setGridLineVisible(false);
            axisTop->setRange(minX, maxX);
            axisTop->setLinePenColor(Qt::black);

            axisRight->setTickAnchor(0);
            axisRight->setTickInterval(deltaY);
            axisRight->setTickType(QValueAxis::TicksDynamic);
            axisRight->setLabelsVisible(false);
            axisRight->setGridLineVisible(false);
            axisRight->setMinorGridLineVisible(false);
            axisRight->setMinorTickCount(4);
            axisRight->setRange(minY, maxY);
            axisRight->setLinePenColor(Qt::black);

            chart.addAxis(axisX, Qt::AlignBottom);
            chart.addAxis(axisY, Qt::AlignLeft);
            chart.addAxis(axisTop, Qt::AlignTop);
            chart.addAxis(axisRight, Qt::AlignRight);

            {
                auto& line = Base::line(VectorType{minX, maxX}, VectorType{0, 0});
                auto pen = line.pen();
                pen.setColor(Qt::black);
                pen.setStyle(Qt::DashLine);
                line.setPen(pen);
                line.attachAxis(axisX);
                line.attachAxis(axisY);

                auto& line1 = Base::line(VectorType{0, 0}, VectorType{minY, maxY});
                line1.setPen(pen);
                line1.attachAxis(axisX);
                line1.attachAxis(axisY);

                pen.setColor(Qt::blue);
                auto& line2 = Base::line(VectorType{0, maxX}, VectorType{0, maxY / numSigma});
                line2.setPen(pen);
                line2.attachAxis(axisX);
                line2.attachAxis(axisY);

                auto& line3 = Base::line(VectorType{minX, 0}, VectorType{maxY / numSigma, 0});
                line3.setPen(pen);
                line3.attachAxis(axisX);
                line3.attachAxis(axisY);
            }
        }
    }
}
