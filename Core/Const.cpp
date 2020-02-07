/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Header/Const.h"
#include "Header/Solve.h"
#include "Header/BasicCalculates.h"
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
    R_MAX = new RealNumber(byte, 10, 9);

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 0;
    ZERO = new RealNumber(byte, 1, 0);

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 1;
    ONE = new RealNumber(byte, 1, 0);

    auto temp = new RealNumber(ONE);
    temp->sign = false;
    MINUS_ONE = temp;

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 2;
    TWO = new RealNumber(byte, 1, 0);

    temp = new RealNumber(TWO);
    temp->sign = false;
    MINUS_TWO = temp;

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 3;
    THREE = new RealNumber(byte, 1, 0);

    temp = new RealNumber(THREE);
    temp->sign = false;
    MINUS_THREE = temp;

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 4;
    FOUR = new RealNumber(byte, 1, 0);

    temp = new RealNumber(FOUR);
    temp->sign = false;
    MINUS_FOUR = temp;
}

Const_1::~Const_1() {
    delete R_MAX;
    delete ZERO;
    delete ONE;
    delete MINUS_ONE;
    delete TWO;
    delete MINUS_TWO;
    delete THREE;
    delete MINUS_THREE;
    delete FOUR;
    delete MINUS_FOUR;
}
/*
 * Consts that need some calculates.
 * Should call new to Const_1 so as to make calculates available.
 */
Const_2::Const_2() {
    //TODO Slow and not accurate
    auto byte = (unsigned char*)malloc(16 * sizeof(char));
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
    PI = new RealNumber(byte, 16, 0);
    E = bisectionMethod(ln_noCheck, const_1->ONE, const_1->TWO, const_1->THREE);

    PI_DIVIDE_TWO = *PI / *const_1->TWO;
    auto temp = new RealNumber(PI_DIVIDE_TWO);
    temp->sign = false;
    MINUS_PI_DIVIDE_TWO = temp;
}

Const_2::~Const_2() {
    delete PI;
    delete E;
    delete PI_DIVIDE_TWO;
    delete MINUS_PI_DIVIDE_TWO;
}