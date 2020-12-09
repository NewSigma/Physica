/*
 * Copyright 2020 WeiBo He.
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
#include <QtGui/QPainter>
#include <QMouseEvent>
#include "Physica/Gui/Plot/DensityPlot.h"
#include "Physica/Logger/LoggerRuntime.h"

namespace Physica::Gui {
    DensityPlot::DensityPlot(int width, int height, QWidget* parent)
            : QWidget(parent)
            , dataImage(QImage(width, height, QImage::Format_Alpha8))
            , colorImage(QImage(width, height, QImage::Format_ARGB32)) {
        setFixedSize(width, height);
        dataImage.fill(Qt::transparent);
        colorImage.fill(Qt::transparent);
        initColorList();
    }

    void DensityPlot::paintEvent(QPaintEvent *event) {
        QPainter painter(this);
        painter.setPen(QPen(Qt::black, 2));
        for(int i = 0; i < width(); i += 50)
            painter.drawLine(i, 0, i, height());
        for(int i=0; i < height(); i += 50)
            painter.drawLine(0, i, width(), i);

        painter.drawImage(0, 0, colorImage);
    }

    void DensityPlot::appendPoint(int x, int y, int radius, unsigned char alpha) {
        if(x < 0 || y < 0 || radius < 0 || x >= width() || y >= height())
            Warning(0, "Encountered invalid argument.");

        QPainter painter(&dataImage);
        painter.setPen(Qt::transparent);

        QRadialGradient gradient(x, y, radius);
        gradient.setColorAt(0,QColor(0, 0, 0, alpha));
        gradient.setColorAt(1,QColor(0, 0, 0, 0));
        painter.setBrush(gradient);
        painter.drawEllipse(QPoint(x, y), radius, radius);
        drawColorImage();
        update();
    }

    void DensityPlot::drawColorImage() {
        for(int row = 0; row < dataImage.height(); row++) {
            const unsigned char* line_data = dataImage.scanLine(row);
            QRgb* colorLine = reinterpret_cast<QRgb*>(colorImage.scanLine(row));
            for(int col = 0; col < dataImage.width(); col++)
                colorLine[col] = colorList[line_data[col]];
        }
    }

    void DensityPlot::initColorList() {
        /*
         * It is suspected that QLinearGradient does not provide a interface to access its color,
         * we can draw it on a image to get the color but a better choice may be add a getter to QLinearGradient.
         */
        QLinearGradient linear(QPoint(0,0),QPoint(255,0));
        linear.setColorAt(0, Qt::blue);
        linear.setColorAt(0.2, Qt::blue);
        linear.setColorAt(0.4, Qt::cyan);
        linear.setColorAt(0.6, Qt::green);
        linear.setColorAt(0.8, Qt::yellow);
        linear.setColorAt(0.95, Qt::red);

        QImage img(256,1,QImage::Format_ARGB32);
        QPainter painter(&img);
        painter.fillRect(img.rect(),linear);

        for(quint32 i = 0;i < 256; i++){
            quint32 alpha = GlobalAlpha / 255.0 * i;
            colorList[i]= (img.pixel(i, 0) & 0x00FFFFFFu) | (alpha << 24U);
        }
    }
}