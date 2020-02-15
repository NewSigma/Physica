#include <iostream>
#include "Core/Header/Const.h"
#include "Core/Header/ElementaryFunction.h"
/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
const Const_1* const_1;
const Const_2* const_2;

int main(int argc, char** argv) {
    const_1 = new Const_1();
    const_2 = new Const_2();
    /*
    QApplication app(argc, argv);
    QApplication::setApplicationName("Physica");
    QApplication::setApplicationVersion("0.0.1");

    QSurfaceFormat format;

    auto main_window = new MainWindow();
    main_window->show();

    int exit_code = QApplication::exec();
    */
    auto a = const_1->getTwo();
    a->power = -1;
    auto result = ln(*a);
    std::cout << *result << std::endl << *const_2->PI;

    delete const_1;
    delete const_2;

    return 0;
}