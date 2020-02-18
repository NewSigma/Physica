/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Header/Const.h"
#include "Header/Solve.h"
#include "Header/ElementaryFunction.h"
#include <malloc.h>

extern const Const_1* const_1;
/*
 * Basic consts that initialize directly.
 */
Const_1::Const_1() {
    GlobalPrecision = 16;

    auto byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 1;
    StepSize = new RealNumber(byte, 1, 1 - GlobalPrecision);

    byte = (unsigned char*)malloc(10 * sizeof(char));
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
    delete StepSize;
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
    PI = bisectionMethod(sin, *const_1->ZERO, *const_1->THREE, *const_1->FOUR);
    E = bisectionMethod(ln_noCheck, *const_1->ONE, *const_1->TWO, *const_1->THREE);

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