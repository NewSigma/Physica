/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtWidgets/QApplication>
#include "Core/Header/Const.h"
#include "Gui/Header/PhysicaMain.h"

const Const_1* const_1;
const Const_2* const_2;

int main(int argc, char** argv) {
    const_1 = new Const_1();
    const_2 = new Const_2();

    QApplication::setApplicationName("Physica");
    QApplication::setApplicationVersion("0.0.1");
    QApplication::setOrganizationName("NewSigma@163.com");

    QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QApplication app(argc, argv);

    new PhysicaMain();

    int exit_code = QApplication::exec();

    delete const_1;
    delete const_2;

    return exit_code;
}