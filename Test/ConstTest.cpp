/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <Const.h>
#include <iostream>
#include "Numerical.h"

void const_test() {
    std::cout << "GlobalPrecision: " << (int) basicConst->GlobalPrecision << '\n';
    std::cout << "ExpectedRelativeError: " << basicConst->getExpectedRelativeError() << '\n';
    std::cout << "StepSize: " << basicConst->getStepSize() << '\n';
    std::cout << "RandomMax: " << basicConst->getR_MAX() << '\n';
    std::cout << "1: " << basicConst->get_1() << '\n';
    std::cout << "-1: " << basicConst->getMinus_1() << '\n';
    std::cout << "2: " << basicConst->get_2() << '\n';
    std::cout << "-2: " << basicConst->getMinus_2() << '\n';
    std::cout << "3: " << basicConst->get_3() << '\n';
    std::cout << "-3: " << basicConst->getMinus_3() << '\n';
    std::cout << "4: " << basicConst->get_4() << '\n';
    std::cout << "-4: " << basicConst->getMinus_4() << '\n';

    std::cout << "Pi: " << mathConst->getPI() << '\n';
    std::cout << "E: " << mathConst->getE() << '\n';
    std::cout << "Pi / 2: " << mathConst->getPI_2() << '\n';
    std::cout << "-Pi / 2: " << mathConst->getMinus_PI_2() << '\n';
}
