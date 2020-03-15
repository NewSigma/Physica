/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#include <QtGui/QGuiApplication>
#include "QScreen"
#include "Header/Settings.h"

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
    setLayout(default_layout);

    list = new QListWidget(this);
    default_layout->addWidget(list);

    themeLabel = new QLabel("Theme", this);
    default_layout->addWidget(themeLabel);

    themeBox = new QComboBox(this);
    default_layout->addWidget(themeBox);

    show();
}