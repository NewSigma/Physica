/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/ComplexNum.h"

ComplexNum::ComplexNum(RealNumber* r, RealNumber* i) {
    real = r;
    imagine = i;
}

ComplexNum::~ComplexNum() {
    delete real;
    delete imagine;
}