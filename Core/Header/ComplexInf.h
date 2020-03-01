/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_COMPLEXINF_H
#define PHYSICA_COMPLEXINF_H

#include "Number.h"

class ComplexInf : public Number {
private:
    bool real_sign;
    bool imagine_sign;
public:
    ComplexInf(bool real_sign, bool imagine_sign);
};

#endif
