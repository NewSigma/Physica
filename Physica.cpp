#include <QtWidgets/QApplication>
#include <iostream>
#include "Core/Header/Const.h"
#include "Gui/Header/MainWindow.h"

/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
const Const_1* const_1;
const Const_2* const_2;

int main(int argc, char** argv) {
    const_1 = new Const_1();
    const_2 = new Const_2();

    QApplication app(argc, argv);
    QApplication::setApplicationName("Physica");
    QApplication::setApplicationVersion("0.0.1");

    new MainWindow();

    int exit_code = QApplication::exec();

    delete const_1;
    delete const_2;

    return exit_code;
}