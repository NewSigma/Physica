/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Physica.h"
#include <QApplication>
#include <iostream>
#include "Const.h"
#include "PhysicaMain.h"

const Const_1* const_1;
const Const_2* const_2;

void init() {
    const_1 = new Const_1();
    const_2 = new Const_2();
}

int main(int argc, char** argv) {
    init();

    QApplication::setApplicationName("Physica");
    QApplication::setApplicationVersion("0.0.1");
    QApplication::setOrganizationName("NewSigma@163.com");

    QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
    QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QApplication app(argc, argv);

    new PhysicaMain();

    int exit_code = QApplication::exec();

    deInit();
    return exit_code;
}

void deInit() {
    delete const_1;
    delete const_2;
}