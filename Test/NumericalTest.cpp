/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/Numerical.h"

void numerical_test() {
    std::cout << "Please input two double numbers. Let the second number be zero to stop the test." << '\n';
    double d1, d2 = 1;
    std::cin >> d1 >> d2;
    while(d2 != 0) {
        Numerical n1(d1);
        Numerical n2(d2);
        auto c1 = n1 + n2;
        std::cout << "d1 + d2 = " << '\n' << c1 << std::endl;
        auto c2 = n1 - n2;
        std::cout << "d1 - d2 = " << '\n' << c2 << std::endl;
        auto c3 = n1 * n2;
        std::cout << "d1 * d2 = " << '\n' << c3 << std::endl;
        auto c4 = n1 / n2;
        std::cout << "d1 / d2 = " << '\n' << c4 << '\n' << "Waiting for next inputs." << std::endl;
        std::cin >> d1 >> d2;
    }
}