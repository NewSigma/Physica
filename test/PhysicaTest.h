/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_PHYSICATEST_H
#define PHYSICA_PHYSICATEST_H

namespace Physica::Test {
    void checkGPU();
    void constTest();
    void elementary_function_test();
    void simpleNet();
    void numericalAddTest(int loop);
    void numericalSubTest(int loop);
    void numericalMulTest(int loop);
    void numericalDivTest(int loop);
    void printElements(const Core::MultiScalar& n);
}

#endif