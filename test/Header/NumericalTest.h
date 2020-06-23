/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_NUMERICALTEST_H
#define PHYSICA_NUMERICALTEST_H

namespace Physica::Core {
    class Scalar;
}

namespace Physica::Test {
    void numericalAddTest(int loop);
    void numericalSubTest(int loop);
    void numericalMulTest(int loop);
    void numericalDivTest(int loop);
    void printElements(const Core::Scalar& n);
}

#endif