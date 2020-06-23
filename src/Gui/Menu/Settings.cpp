/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#include <QtGui/QGuiApplication>
#include <QScreen>
#include "Physica/Gui/Settings.h"

Settings::Settings(PhysicaMain* parent) : parent(parent) {
    //Basic settings
    /* Resize */ {
        QRect primaryScreenRec = QGuiApplication::primaryScreen()->geometry();
        resize(primaryScreenRec.width() * 8 / 10, primaryScreenRec.height() * 8 / 10);
    }
    setAttribute(Qt::WA_DeleteOnClose);
    setWindowTitle("Settings");
    //Main
    default_layout = new QHBoxLayout(this);

    list = new QListWidget(this);
    list->setFixedSize(width() / 5, height());
    default_layout->addWidget(list);

    page = new SettingsPage(this);
    default_layout->addWidget(page);

    show();
}