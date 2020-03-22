/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Header/Const.h"
#include "Header/Solve.h"
#include "Header/ElementaryFunction.h"
#include <malloc.h>
/*
 * Basic consts that initialize directly.
 */
Const_1::Const_1() {
    GlobalPrecision = 16;

    auto byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 1;
    //StepSize is also is the smallest relative error that is acceptable.
    stepSize = new Numerical(byte, 1, 1 - GlobalPrecision);

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
    R_MAX = new Numerical(byte, 10, 9);

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 0;
    _0 = new Numerical(byte, 1, 0);

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 1;
    _1 = new Numerical(byte, 1, 0);

    auto temp = new Numerical(_1);
    temp->sign = false;
    Minus_1 = temp;

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 2;
    _2 = new Numerical(byte, 1, 0);

    temp = new Numerical(_2);
    temp->sign = false;
    Minus_2 = temp;

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 3;
    _3 = new Numerical(byte, 1, 0);

    temp = new Numerical(_3);
    temp->sign = false;
    Minus_3 = temp;

    byte = (unsigned char*)malloc(sizeof(char));
    byte[0] = 4;
    _4 = new Numerical(byte, 1, 0);

    temp = new Numerical(_4);
    temp->sign = false;
    Minus_4 = temp;
}

Const_1::~Const_1() {
    delete stepSize;
    delete R_MAX;
    delete _0;
    delete _1;
    delete Minus_1;
    delete _2;
    delete Minus_2;
    delete _3;
    delete Minus_3;
    delete _4;
    delete Minus_4;
}
/*
 * Consts that need some calculates.
 * Should call new to Const_1 so as to make calculates available.
 */
Const_2::Const_2() {
    stepSize = new RealNum(new Numerical(const_1->stepSize));
    _0 = new RealNum(getZero());
    _1 = new RealNum(getOne());
    _2 = new RealNum(getTwo());
    //TODO Slow and not accurate
    PI = new Numerical(bisectionMethod(sin, *const_1->_0, *const_1->_3, *const_1->_4));
    E = exp(*const_1->_1);

    PI_2 = *PI / *const_1->_2;
    auto temp = new Numerical(PI_2);
    temp->sign = false;
    Minus_PI_2 = temp;
}

Const_2::~Const_2() {
    delete stepSize;
    delete _0;
    delete _1;
    delete _2;
    delete PI;
    delete E;
    delete PI_2;
    delete Minus_PI_2;
}