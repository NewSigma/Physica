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
#include <QtGui/QGuiApplication>
#include <QScreen>
#include "Physica/Gui/Settings.h"

namespace Physica::Gui {
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
}
