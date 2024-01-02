/*
 * Copyright 2021 WeiBo He.
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
#include <QtGui/QGuiApplication>
#include "Physica/Gui/Plot/Plot3D.h"

namespace Physica::Gui {
    Plot3D::Plot3D(QWidget* parent)
            : QWidget(parent)
            , vLayout(new QVBoxLayout())
            , surface(new Q3DSurface()) {
        if (parent == nullptr) {
            QRect primaryScreenRec = QGuiApplication::primaryScreen()->geometry();
            resize(primaryScreenRec.width() / 2, primaryScreenRec.height() / 1.6);
        }
        setAttribute(Qt::WA_DeleteOnClose);
        setLayout(vLayout);

        surface->activeTheme()->setType(Q3DTheme::ThemePrimaryColors);
        surface->axisX()->setLabelAutoRotation(30);
        surface->axisY()->setLabelAutoRotation(90);
        surface->axisZ()->setLabelAutoRotation(30);

        auto* widget = QWidget::createWindowContainer(surface, this);
        widget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        layout()->addWidget(widget);
    }

    QLinearGradient Plot3D::makeDefaultGrad() {
        QLinearGradient gr;
        gr.setColorAt(0.0, Qt::blue);
        gr.setColorAt(0.35, Qt::cyan);
        gr.setColorAt(0.5, Qt::green);
        gr.setColorAt(0.65, Qt::yellow);
        gr.setColorAt(1.0, Qt::red);
        return gr;
    }
}
