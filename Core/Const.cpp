/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Header/Const.h"

extern const Const_1* const_1;
/*
 * Basic consts that initialize directly.
 */
Const_1::Const_1() {
    MachinePrecision = 16;
    R_MAX = new RealNumberA(new unsigned char[10]{2,1,4,7,4,8,3,6,4,7}, 10, 9, true);
    ZERO = new RealNumberA(new unsigned char[1]{0}, 1, 0, true);
    ONE = new RealNumberA(new unsigned char[1]{1}, 1, 0, true);
    MINUS_ONE = new RealNumberA(new unsigned char[1]{1}, 1, 0, false);
    TWO = new RealNumberA(new unsigned char[1]{2}, 1, 0, true);
    MINUS_TWO = new RealNumberA(new unsigned char[1]{2}, 1, 0, false);
    PI = new RealNumberA(new unsigned char[16]{3,1,4,1,5,9,2,6,5,3,5,8,9,7,9,3}, 16, 0, true);
}

Const_1::~Const_1() {
    delete R_MAX;
    delete ZERO;
    delete ONE;
    delete MINUS_ONE;
    delete TWO;
    delete MINUS_TWO;
    delete PI;
}
/*
 * Consts that need some calculates.
 * Should call new to Const_1 so as to make calculates available.
 */
Const_2::Const_2() {
    PI_DIVIDE_TWO = *const_1->PI / *const_1->TWO;
    auto temp = new RealNumberA(PI_DIVIDE_TWO);
    temp->sign = false;
    MINUS_PI_DIVIDE_TWO = temp;
}

Const_2::~Const_2() {
    delete PI_DIVIDE_TWO;
    delete MINUS_PI_DIVIDE_TWO;
}