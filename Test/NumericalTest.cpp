/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/Numerical.h"
#include <random>
#include <iomanip>

using namespace Physica::Core;

namespace Physica::Test {
    void numericalAddTest(int loop) {
        double d;
        std::default_random_engine engine(clock());
        std::cout << "Performing add test" << '\n' << std::setprecision(10);
        for(int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            Numerical a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            Numerical b(d_b);

            double_extract expect{d_a + d_b};
            double_extract result{double(a + b)};
            if(expect.value != result.value) {
                std::cout << "Performing add test " << d_a << " + " << d_b << '\n';
                std::cout << "Performing add test " << "--Failed\n"
                        << "low:\t" << expect.structure.low << '\t' << result.structure.low
                        << "\nhigh:\t" << expect.structure.high << '\t' << result.structure.high
                        << "\nexp:\t" << expect.structure.exp << '\t' << result.structure.exp
                        << "\nsign:\t" << expect.structure.sign << '\t' << result.structure.sign;
                return;
            }
        }
        std::cout << "Performing add test: " << " --Passed" << '\n' << std::setprecision(6);
    }

    void numericalSubTest(int loop) {
        double d;
        std::default_random_engine engine(clock());
        std::cout << "Performing add test" << '\n' << std::setprecision(10);
        for(int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            Numerical a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            Numerical b(d_b);

            double_extract expect{d_a - d_b};
            double_extract result{double(a - b)};
            if(expect.value != result.value) {
                std::cout << "Performing add test " << d_a << " - " << d_b << '\n';
                std::cout << "Performing add test " << "--Failed\n"
                          << "low:\t" << expect.structure.low << '\t' << result.structure.low
                          << "\nhigh:\t" << expect.structure.high << '\t' << result.structure.high
                          << "\nexp:\t" << expect.structure.exp << '\t' << result.structure.exp
                          << "\nsign:\t" << expect.structure.sign << '\t' << result.structure.sign;
                return;
            }
        }
        std::cout << "Performing add test: " << " --Passed" << '\n' << std::setprecision(6);
    }

    void numericalMulTest(int loop) {
        double d;
        std::default_random_engine engine(clock());
        std::cout << "Performing add test" << '\n' << std::setprecision(10);
        for(int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            Numerical a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            Numerical b(d_b);

            double_extract expect{d_a * d_b};
            double_extract result{double(a * b)};
            if(expect.value != result.value) {
                std::cout << "Performing add test " << d_a << " * " << d_b << '\n';
                std::cout << "Performing add test " << "--Failed\n"
                          << "low:\t" << expect.structure.low << '\t' << result.structure.low
                          << "\nhigh:\t" << expect.structure.high << '\t' << result.structure.high
                          << "\nexp:\t" << expect.structure.exp << '\t' << result.structure.exp
                          << "\nsign:\t" << expect.structure.sign << '\t' << result.structure.sign;
                return;
            }
        }
        std::cout << "Performing add test: " << " --Passed" << '\n' << std::setprecision(6);
    }

    void numericalDivTest(int loop) {
        double d;
        std::default_random_engine engine(clock());
        std::cout << "Performing add test" << '\n' << std::setprecision(10);
        for(int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            Numerical a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            while(d_b == 0) {
                d = 1 - 1.4 * d * d;
                d_b = d * engine();
            }
            Numerical b(d_b);

            double_extract expect{d_a / d_b};
            double_extract result{double(a / b)};
            if(expect.value != result.value) {
                std::cout << "Performing add test " << d_a << " / " << d_b << '\n';
                std::cout << "Performing add test " << "--Failed\n"
                          << "low:\t" << expect.structure.low << '\t' << result.structure.low
                          << "high:\t" << expect.structure.high << '\t' << result.structure.high
                          << "exp:\t" << expect.structure.exp << '\t' << result.structure.exp
                          << "sign:\t" << expect.structure.sign << '\t' << result.structure.sign;
                return;
            }
        }
        std::cout << "Performing add test: " << " --Passed" << '\n' << std::setprecision(6);
    }

    void printElements(const Numerical& n) {
        int size = n.getSize();
        for(int i = 0; i < size; ++i)
            std::cout << n[i] << ' ';
    }
}