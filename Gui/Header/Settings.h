/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SETTINGS_H
#define PHYSICA_SETTINGS_H

#include "PhysicaMain.h"
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QLabel>
#include <QtWidgets/QComboBox>
#include "QDialog"

class Settings : public QDialog {
    Q_OBJECT

    PhysicaMain* parent;
    QHBoxLayout* default_layout;
    QListWidget* list;
    QLabel* themeLabel;
    QComboBox* themeBox;
public:
    explicit Settings(PhysicaMain* parent);
};


#endif
