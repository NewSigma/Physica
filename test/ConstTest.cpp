/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <iostream>
#include <Physica/Core/MultiPrecition/Scalar.h>
#include "PhysicaTest.h"

using namespace Physica::Core;

namespace Physica::Test {
    void constTest() {
        std::cout << "GlobalPrecision: " << (int) GlobalPrecision << '\n';
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

        /* test Pi */ {
            auto Pi = MathConst::getInstance().getPI();
            std::cout << "Testing Pi" << '\n';
            if(Pi.getSize() == 5 && Pi[0] == 11424456171093639400UL && Pi[1] == 11820040416388919749UL
                && Pi[2] == 1376283091369227076 && Pi[3] == 2611923443488327891 && Pi[4] == 3)
                std::cout << "Testing Pi --Passed" << '\n';
            else {
                std::cout << "Testing Pi --Failed" << '\n' << "Pi: ";
                printElements(MathConst::getInstance().getPI());
            }
        }
        /* test E */ {
            auto E = MathConst::getInstance().getE();
            std::cout << "Testing E" << '\n';
            if(E.getSize() == 4 && E[0] == 7126689189968796226
               && E[1] == 13794904443024896967UL && E[2] == 13249961062380153450UL && E[3] == 2)
                std::cout << "Testing E --Passed" << '\n';
            else {
                std::cout << "Testing E --Failed" << '\n' << "Pi: ";
                printElements(MathConst::getInstance().getE());
            }
        }
        std::cout << "Pi / 2: " << MathConst::getInstance().getPI_2() << '\n';
        std::cout << "-Pi / 2: " << MathConst::getInstance().getMinus_PI_2() << '\n';
    }
}