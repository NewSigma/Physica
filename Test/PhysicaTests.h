/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PHYSICATESTS_H
#define PHYSICA_PHYSICATESTS_H

namespace Physica::Core {
    class Numerical;
}

using Physica::Core::Numerical;

namespace Physica::Test {
    void checkGPU();
    void constTest();
    void elementary_function_test();
    void numericalTest();
    void printElements(const Numerical& n);
    void simpleNet();
}

#endif