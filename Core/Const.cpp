/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Header/Const.h"
#include <malloc.h>

extern const Const_1* const_1;
/*
 * Basic consts that initialize directly.
 */
Const_1::Const_1() {
    MachinePrecision = 16;
    auto byte = (unsigned char*)malloc(10 * sizeof(char));
    byte[0] = 2;
    byte[1] = 1;
    byte[2] = 4;
    byte[3] = 7;
    byte[4] = 4;
    byte[5] = 8;
    byte[6] = 3;
    byte[7] = 6;
    byte[8] = 4;
    byte[9] = 7;
    R_MAX = new RealNumberA(byte, 10, 9, true);

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 0;
    ZERO = new RealNumberA(byte, 1, 0, true);

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 1;
    RealNumberA* temp;
    ONE = new RealNumberA(byte, 1, 0, true);

    temp = new RealNumberA(ONE);
    temp->sign = false;
    MINUS_ONE = temp;

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 2;
    TWO = new RealNumberA(byte, 1, 0, true);

    temp = new RealNumberA(TWO);
    temp->sign = false;
    MINUS_TWO = temp;

    byte = (unsigned char*)malloc(16 * sizeof(char));
    byte[0] = 3;
    byte[1] = 1;
    byte[2] = 4;
    byte[3] = 1;
    byte[4] = 5;
    byte[5] = 9;
    byte[6] = 2;
    byte[7] = 6;
    byte[8] = 5;
    byte[9] = 3;
    byte[10] = 5;
    byte[11] = 8;
    byte[12] = 9;
    byte[13] = 7;
    byte[14] = 9;
    byte[15] = 3;
    PI = new RealNumberA(byte, 16, 0, true);
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