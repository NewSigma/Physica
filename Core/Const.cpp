/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Const.h"
#include "Numerical.h"
#include "climits"
#include "RealNum.h"
/*
 * Basic consts that initialize directly.
 */
BasicConst::BasicConst() : GlobalPrecision(4) {
    auto byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
    byte[0] = 100;
    plotPoints = new Numerical(byte, 1, 0);

    byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
    byte[0] = 1;
    expectedRelativeError = new Numerical(byte, 1, 1 - GlobalPrecision);

    byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
    byte[0] = 1;
    //Value (- GlobalPrecision / 2) still need a proof.
    stepSize = new Numerical(byte, 1, - GlobalPrecision / 2);

    R_MAX = new Numerical(2147483647);

    byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
    byte[0] = 0;
    _0 = new Numerical(byte, 1, 0);

    byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
    byte[0] = 1;
    _1 = new Numerical(byte, 1, 0);
    Minus_1 = new Numerical(-*_1);

    byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
    byte[0] = 2;
    _2 = new Numerical(byte, 1, 0);
    Minus_2 = new Numerical(-*_2);

    byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
    byte[0] = 3;
    _3 = new Numerical(byte, 1, 0);
    Minus_3 = new Numerical(-*_3);

    byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
    byte[0] = 4;
    _4 = new Numerical(byte, 1, 0);
    Minus_4 = new Numerical(-*_4);
}

BasicConst::~BasicConst() {
    delete expectedRelativeError;
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
MathConst::MathConst() {
    stepSize = new RealNum(basicConst->getStepSize());
    _0 = new RealNum(getZero());
    _1 = new RealNum(getOne());
    _2 = new RealNum(getTwo());
    //0.31 is the big approximation of ln(2) / ln(10)
    PI = new Numerical(calcPI((int)(LONG_WIDTH * basicConst->getGlobalPrecision() * 0.31) + 1));
    E = new Numerical(exp(basicConst->get_1()));

    PI_2 = new Numerical(*PI / basicConst->get_2());
    Minus_PI_2 = new Numerical(-*PI_2);
}

MathConst::~MathConst() {
    delete stepSize;
    delete _0;
    delete _1;
    delete _2;
    delete PI;
    delete E;
    delete PI_2;
    delete Minus_PI_2;
}
/*
 * precision is the number of effective digits in decimal.
 * Reference:
 * http://www.pi314.net/eng/salamin.php
 * https://blog.csdn.net/liangbch/article/details/78724041
 */
Numerical MathConst::calcPI(int precision) {
    Numerical a = getOne();
    Numerical x = getOne();
    Numerical b = reciprocal(sqrt(basicConst->get_2()));
    Numerical c = reciprocal(basicConst->get_4());

    int goal = 1;
    while(goal < precision) {
        Numerical y(a);
        a = (a + b) / basicConst->get_2();
        b = sqrt(b * y);
        y -= a;
        c -= y * y * x;
        x *= basicConst->get_2();
        goal *= 2;
    }
    a = (a + b) / basicConst->get_2();
    return a * a / c;
}