/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/ComplexNum.h"

ComplexNum::ComplexNum(class RealNumber* r, class RealNumber* i) {
    type = Number::ComplexNumber;
    real = r;
    imagine = i;
}

ComplexNum::~ComplexNum() {
    delete real;
    delete imagine;
}