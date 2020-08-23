/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

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
