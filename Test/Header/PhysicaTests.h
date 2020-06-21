/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PHYSICATESTS_H
#define PHYSICA_PHYSICATESTS_H

namespace Physica::Core {
    class Scalar;
}

using Physica::Core::Scalar;

namespace Physica::Test {
    void checkGPU();
    void constTest();
    void elementary_function_test();
    void simpleNet();
}

#endif