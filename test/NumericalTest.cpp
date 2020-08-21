/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <random>
#include <iomanip>
#include <iostream>
#include <Physica/Physica.h>
#include <Physica/Core/MultiPrecition/Scalar.h>

using Physica::Core::MultiScalar;

void numericalAddTest(int loop) {
    double d;
    std::default_random_engine engine(clock());
    std::cout << "Performing add test" << '\n' << std::setprecision(10);
    for(int i = 0; i < loop; ++i) {
        d = 1 - 1.4 * d * d;
        double d_a = d * engine();
        MultiScalar a(d_a);

        d = 1 - 1.4 * d * d;
        double d_b = d * engine();
        MultiScalar b(d_b);

        double_extract expect{d_a + d_b};
        double_extract result{double(a + b)};
        if((expect.sign != result.sign) //NOLINT
           || (expect.exp != result.exp) //NOLINT
           || (expect.high != result.high) //NOLINT
           || abs(static_cast<int>(expect.low) - static_cast<int>(result.low)) > 1) {
            std::cout << "Performing add test " << d_a << " + " << d_b << '\n';
            std::cout << "Performing add test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                      << "low:\t" << expect.low << '\t' << result.low << '\n'
                      << "\nhigh:\t" << expect.high << '\t' << result.high << '\n'
                      << "\nexp:\t" << expect.exp << '\t' << result.exp << '\n'
                      << "\nsign:\t" << expect.sign << '\t' << result.sign << '\n';
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
        MultiScalar a(d_a);

        d = 1 - 1.4 * d * d;
        double d_b = d * engine();
        MultiScalar b(d_b);

        double_extract expect{d_a - d_b};
        double_extract result{double(a - b)};
        if((expect.sign != result.sign) //NOLINT
           || (expect.exp != result.exp) //NOLINT
           || (expect.high != result.high) //NOLINT
           || abs(static_cast<int>(expect.low) - static_cast<int>(result.low)) > 1) {
            std::cout << "Performing sub test " << d_a << " - " << d_b << '\n';
            std::cout << "Performing sub test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                      << "low:\t" << expect.low << '\t' << result.low << '\n'
                      << "\nhigh:\t" << expect.high << '\t' << result.high << '\n'
                      << "\nexp:\t" << expect.exp << '\t' << result.exp << '\n'
                      << "\nsign:\t" << expect.sign << '\t' << result.sign << '\n';
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
        MultiScalar a(d_a);

        d = 1 - 1.4 * d * d;
        double d_b = d * engine();
        MultiScalar b(d_b);

        double_extract expect{d_a * d_b};
        double_extract result{double(a * b)};
        if((expect.sign != result.sign) //NOLINT
           || (expect.exp != result.exp) //NOLINT
           || (expect.high != result.high) //NOLINT
           || abs(static_cast<int>(expect.low) - static_cast<int>(result.low)) > 1) {
            std::cout << "Performing mul test " << d_a << " * " << d_b << '\n';
            std::cout << "Performing mul test " "--Failed (" << (i + 1) << '/' << loop << ")\n"
                      << "low:\t" << expect.low << '\t' << result.low << '\n'
                      << "\nhigh:\t" << expect.high << '\t' << result.high << '\n'
                      << "\nexp:\t" << expect.exp << '\t' << result.exp << '\n'
                      << "\nsign:\t" << expect.sign << '\t' << result.sign << '\n';
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
        MultiScalar a(d_a);

        d = 1 - 1.4 * d * d;
        double d_b = d * engine();
        while(d_b == 0) {
            d = 1 - 1.4 * d * d;
            d_b = d * engine();
        }
        MultiScalar b(d_b);

        double_extract expect{d_a / d_b};
        double_extract result{double(a / b)};
        if((expect.sign != result.sign) //NOLINT
           || (expect.exp != result.exp) //NOLINT
           || (expect.high != result.high) //NOLINT
           || abs(static_cast<int>(expect.low) - static_cast<int>(result.low)) > 1) {
            std::cout << "Performing div test " << d_a << " / " << d_b << '\n';
            std::cout << "Performing div test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                      << "low:\t" << expect.low << '\t' << result.low << '\n'
                      << "high:\t" << expect.high << '\t' << result.high << '\n'
                      << "exp:\t" << expect.exp << '\t' << result.exp << '\n'
                      << "sign:\t" << expect.sign << '\t' << result.sign << '\n';
            return;
        }
    }
    std::cout << "Performing div test: " << " --Passed" << '\n' << std::setprecision(6);
}

int main(int argc, char** argv) {
    initPhysica();
    numericalAddTest(50);
    numericalSubTest(50);
    numericalMulTest(50);
    numericalDivTest(50);
    deInitPhysica();
    return 0;
}