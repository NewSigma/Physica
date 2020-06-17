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
            if((expect.structure.sign != result.structure.sign) //NOLINT
               || (expect.structure.exp != result.structure.exp) //NOLINT
               || (expect.structure.high != result.structure.high) //NOLINT
               || abs(static_cast<int>(expect.structure.low) - static_cast<int>(result.structure.low)) > 1) {
                std::cout << "Performing add test " << d_a << " + " << d_b << '\n';
                std::cout << "Performing add test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                        << "low:\t" << expect.structure.low << '\t' << result.structure.low << '\n'
                        << "\nhigh:\t" << expect.structure.high << '\t' << result.structure.high << '\n'
                        << "\nexp:\t" << expect.structure.exp << '\t' << result.structure.exp << '\n'
                        << "\nsign:\t" << expect.structure.sign << '\t' << result.structure.sign << '\n';
                return;
            }
        }
        std::cout << "Performing add test: " << " --Passed" << '\n' << std::setprecision(6);
    }

    void numericalSubTest(int loop) {
        double d;
        std::default_random_engine engine(clock());
        std::cout << "Performing sub test" << '\n' << std::setprecision(10);
        for(int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            Numerical a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            Numerical b(d_b);

            double_extract expect{d_a - d_b};
            double_extract result{double(a - b)};
            if((expect.structure.sign != result.structure.sign) //NOLINT
               || (expect.structure.exp != result.structure.exp) //NOLINT
               || (expect.structure.high != result.structure.high) //NOLINT
               || abs(static_cast<int>(expect.structure.low) - static_cast<int>(result.structure.low)) > 1) {
                std::cout << "Performing sub test " << d_a << " - " << d_b << '\n';
                std::cout << "Performing sub test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                          << "low:\t" << expect.structure.low << '\t' << result.structure.low << '\n'
                          << "\nhigh:\t" << expect.structure.high << '\t' << result.structure.high << '\n'
                          << "\nexp:\t" << expect.structure.exp << '\t' << result.structure.exp << '\n'
                          << "\nsign:\t" << expect.structure.sign << '\t' << result.structure.sign << '\n';
                return;
            }
        }
        std::cout << "Performing sub test: " << " --Passed" << '\n' << std::setprecision(6);
    }

    void numericalMulTest(int loop) {
        double d;
        std::default_random_engine engine(clock());
        std::cout << "Performing mul test" << '\n' << std::setprecision(10);
        for(int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            Numerical a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            Numerical b(d_b);

            double_extract expect{d_a * d_b};
            double_extract result{double(a * b)};
            if((expect.structure.sign != result.structure.sign) //NOLINT
               || (expect.structure.exp != result.structure.exp) //NOLINT
               || (expect.structure.high != result.structure.high) //NOLINT
               || abs(static_cast<int>(expect.structure.low) - static_cast<int>(result.structure.low)) > 1) {
                std::cout << "Performing mul test " << d_a << " * " << d_b << '\n';
                std::cout << "Performing mul test " "--Failed (" << (i + 1) << '/' << loop << ")\n"
                          << "low:\t" << expect.structure.low << '\t' << result.structure.low << '\n'
                          << "\nhigh:\t" << expect.structure.high << '\t' << result.structure.high << '\n'
                          << "\nexp:\t" << expect.structure.exp << '\t' << result.structure.exp << '\n'
                          << "\nsign:\t" << expect.structure.sign << '\t' << result.structure.sign << '\n';
                return;
            }
        }
        std::cout << "Performing mul test: " << " --Passed" << '\n' << std::setprecision(6);
    }

    void numericalDivTest(int loop) {
        double d;
        std::default_random_engine engine(clock());
        std::cout << "Performing div test" << '\n' << std::setprecision(10);
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
            if((expect.structure.sign != result.structure.sign) //NOLINT
               || (expect.structure.exp != result.structure.exp) //NOLINT
               || (expect.structure.high != result.structure.high) //NOLINT
               || abs(static_cast<int>(expect.structure.low) - static_cast<int>(result.structure.low)) > 1) {
                std::cout << "Performing div test " << d_a << " / " << d_b << '\n';
                std::cout << "Performing div test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                          << "low:\t" << expect.structure.low << '\t' << result.structure.low << '\n'
                          << "high:\t" << expect.structure.high << '\t' << result.structure.high << '\n'
                          << "exp:\t" << expect.structure.exp << '\t' << result.structure.exp << '\n'
                          << "sign:\t" << expect.structure.sign << '\t' << result.structure.sign << '\n';
                return;
            }
        }
        std::cout << "Performing div test: " << " --Passed" << '\n' << std::setprecision(6);
    }

    void printElements(const Numerical& n) {
        int size = n.getSize();
        for(int i = 0; i < size; ++i)
            std::cout << n[i] << ' ';
    }
}