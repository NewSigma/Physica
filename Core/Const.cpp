/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Const.h"
#include "Numerical.h"
#include "RealNum.h"

namespace Physica::Core {
    BasicConst* BasicConst::instance = nullptr;
    MathConst* MathConst::instance = nullptr;
    /*
     * Basic consts that initialize directly.
     */
    BasicConst::BasicConst() : GlobalPrecision(4), MaxPower(16) {
        plotPoints = new Numerical(static_cast<SignedNumericalUnit>(20));
        auto byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
        byte[0] = 1;
        expectedRelativeError = new Numerical(byte, 1, 1 - GlobalPrecision);
        byte = reinterpret_cast<unsigned long*>(malloc(sizeof(long)));
        byte[0] = 1;
        //Value (- GlobalPrecision / 2) still need a proof.
        stepSize = new Numerical(byte, 1, - GlobalPrecision / 2);
        R_MAX = new Numerical(static_cast<SignedNumericalUnit>(2147483647));
        _0 = new Numerical(static_cast<SignedNumericalUnit>(0));
        _1 = new Numerical(static_cast<SignedNumericalUnit>(1));
        Minus_1 = new Numerical(static_cast<SignedNumericalUnit>(-1));
        _2 = new Numerical(static_cast<SignedNumericalUnit>(2));
        Minus_2 = new Numerical(static_cast<SignedNumericalUnit>(-2));
        _3 = new Numerical(static_cast<SignedNumericalUnit>(3));
        Minus_3 = new Numerical(static_cast<SignedNumericalUnit>(-3));
        _4 = new Numerical(static_cast<SignedNumericalUnit>(4));
        Minus_4 = new Numerical(static_cast<SignedNumericalUnit>(-4));
        _10 = new Numerical(static_cast<SignedNumericalUnit>(10));
    }

    BasicConst::~BasicConst() {
        delete plotPoints;
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
        delete _10;
    }

    void BasicConst::init() {
        if(instance == nullptr)
            instance = new BasicConst();
    }

    void BasicConst::deInit() {
        if(instance != nullptr) {
            delete instance;
            instance = nullptr;
        }
    }
    /*
     * Consts that need some calculates.
     * Should call new to Const_1 so as to make calculates available.
     */
    MathConst::MathConst() {
        stepSize = new RealNum(BasicConst::getInstance().getStepSize());
        _0 = new RealNum(getZero());
        _1 = new RealNum(getOne());
        _2 = new RealNum(getTwo());
        //0.31 is the big approximation of ln(2) / ln(10)
        PI = new Numerical(calcPI(
                static_cast<int>(static_cast<double>(NumericalUnitWidth) * BasicConst::getInstance().GlobalPrecision * 0.31) + 1));
        E = new Numerical(exp(BasicConst::getInstance().get_1()));

        PI_2 = new Numerical(*PI / BasicConst::getInstance().get_2());
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

    void MathConst::init() {
        if(instance == nullptr)
            instance = new MathConst();
    }

    void MathConst::deInit() {
        if(instance != nullptr) {
            delete instance;
            instance = nullptr;
        }
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
        Numerical b = reciprocal(sqrt(BasicConst::getInstance().get_2()));
        Numerical c = reciprocal(BasicConst::getInstance().get_4());

        int goal = 1;
        while(goal < precision) {
            Numerical y(a);
            a = (a + b) >> 1U;
            b = sqrt(b * y);
            y -= a;
            c -= y * y * x;
            x *= BasicConst::getInstance().get_2();
            goal *= 2;
        }
        a = (a + b) >> 1U;
        return a * a / c;
    }
}
