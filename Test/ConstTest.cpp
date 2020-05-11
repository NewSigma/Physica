/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <Core/Header/Const.h>
#include <iostream>
#include "Core/Header/Numerical.h"

using namespace Physica::Core;

void constTest() {
    std::cout << "GlobalPrecision: " << (int) BasicConst::getInstance().GlobalPrecision << '\n';
    std::cout << "ExpectedRelativeError: " << BasicConst::getInstance().getExpectedRelativeError() << '\n';
    std::cout << "StepSize: " << BasicConst::getInstance().getStepSize() << '\n';
    std::cout << "RandomMax: " << BasicConst::getInstance().getR_MAX() << '\n';
    std::cout << "1: " << BasicConst::getInstance().get_1() << '\n';
    std::cout << "-1: " << BasicConst::getInstance().getMinus_1() << '\n';
    std::cout << "2: " << BasicConst::getInstance().get_2() << '\n';
    std::cout << "-2: " << BasicConst::getInstance().getMinus_2() << '\n';
    std::cout << "3: " << BasicConst::getInstance().get_3() << '\n';
    std::cout << "-3: " << BasicConst::getInstance().getMinus_3() << '\n';
    std::cout << "4: " << BasicConst::getInstance().get_4() << '\n';
    std::cout << "-4: " << BasicConst::getInstance().getMinus_4() << '\n';

    std::cout << "Pi: " << MathConst::getInstance().getPI() << '\n';
    std::cout << "E: " << MathConst::getInstance().getE() << '\n';
    std::cout << "Pi / 2: " << MathConst::getInstance().getPI_2() << '\n';
    std::cout << "-Pi / 2: " << MathConst::getInstance().getMinus_PI_2() << '\n';
}
