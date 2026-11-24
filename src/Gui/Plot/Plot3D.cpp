/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include <QtGui/QGuiApplication>
#include "Physica/Gui/Plot/Plot3D.h"

using namespace Physica;

Plot3D::Plot3D(QWidget* parent)
        : QWidget(parent)
        , vLayout(new QVBoxLayout())
        , surface(new Q3DSurface()) {
    if (parent == nullptr) {
        QRect primaryScreenRec = QGuiApplication::primaryScreen()->geometry();
        resize(primaryScreenRec.width() / 2, int(double(primaryScreenRec.height()) / 1.6));
    }
    setAttribute(Qt::WA_DeleteOnClose);
    setLayout(vLayout);

    surface->activeTheme()->setType(Q3DTheme::ThemePrimaryColors);
    getAxisX()->setLabelAutoRotation(30);
    getAxisY()->setLabelAutoRotation(90);
    getAxisZ()->setLabelAutoRotation(30);

    getAxisX()->setTitleVisible(true);
    getAxisY()->setTitleVisible(true);
    getAxisZ()->setTitleVisible(true);

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
