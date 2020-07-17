/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_SETTINGS_H
#define PHYSICA_SETTINGS_H

#include "PhysicaMain.h"
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QListWidget>
#include "QDialog"
#include "SettingsPage.h"

namespace Physica::Gui {
    class Settings : public QDialog {
        Q_OBJECT

        PhysicaMain* parent;
        QHBoxLayout* default_layout;
        QListWidget* list;
        SettingsPage* page;
    public:
        explicit Settings(PhysicaMain* parent);
    };

}

#endif
