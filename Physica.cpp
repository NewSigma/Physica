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
#include <QApplication>
#include "Physica/PhysicaInit.h"
#include "Physica/Gui/PhysicaMain.h"

int main(int argc, char** argv) {
    initPhysica();

    QApplication::setApplicationName("Physica");
    QApplication::setApplicationVersion("0.0.1");
    QApplication::setOrganizationName("NewSigma@163.com");

    QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QApplication app(argc, argv);

    new Physica::Gui::PhysicaMain();

    int exit_code = QApplication::exec();

    deInitPhysica();
    return exit_code;
}