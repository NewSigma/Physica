/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../../Header/ComplexInf.h"

ComplexInf::ComplexInf(bool rsign, bool isign) {
    type = Number::ComplexInf;
    real_sign = rsign;
    imagine_sign = isign;
}