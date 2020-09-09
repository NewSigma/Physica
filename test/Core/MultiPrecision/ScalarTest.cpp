/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <random>
#include <iostream>
#include <Physica/PhysicaInit.h>
#include <Physica/Core/MultiPrecition/Scalar.h>

using Physica::Core::MultiScalar;
constexpr unsigned int iterateCount = 50;
static std::default_random_engine engine(clock());
/*!
 * Test operator+(), return true if passed.
 */
bool numericalAddTest(unsigned int loop) {
    double d;
    for(unsigned int i = 0; i < loop; ++i) {
        d = 1 - 1.4 * d * d;
        double d_a = d * engine();
        MultiScalar a(d_a);

        d = 1 - 1.4 * d * d;
        double d_b = d * engine();
        MultiScalar b(d_b);

        double_extract expect{d_a + d_b};
        double_extract result{double(a + b)};
        if((expect.sign != static_cast<unsigned int>(result.sign))
           || (expect.exp != static_cast<unsigned int>(result.exp))
           || (expect.high != static_cast<unsigned int>(result.high))
           || abs(static_cast<int>(expect.low) - static_cast<int>(result.low)) > 1) {
            std::cout << "Performing add test " << d_a << " + " << d_b << '\n';
            std::cout << "Performing add test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                      << "low:\t" << expect.low << '\t' << result.low << '\n'
                      << "\nhigh:\t" << expect.high << '\t' << result.high << '\n'
                      << "\nexp:\t" << expect.exp << '\t' << result.exp << '\n'
                      << "\nsign:\t" << expect.sign << '\t' << result.sign << '\n';
            return false;
        }
    }
    return true;
}
/*!
 * Test operator-(), return true if passed.
 */
bool numericalSubTest(unsigned int loop) {
    double d;
    for(unsigned int i = 0; i < loop; ++i) {
        d = 1 - 1.4 * d * d;
        double d_a = d * engine();
        MultiScalar a(d_a);

        d = 1 - 1.4 * d * d;
        double d_b = d * engine();
        MultiScalar b(d_b);

        double_extract expect{d_a - d_b};
        double_extract result{double(a - b)};
        if((expect.sign != static_cast<unsigned int>(result.sign))
           || (expect.exp != static_cast<unsigned int>(result.exp))
           || (expect.high != static_cast<unsigned int>(result.high))
           || abs(static_cast<int>(expect.low) - static_cast<int>(result.low)) > 1) {
            std::cout << "Performing sub test " << d_a << " - " << d_b << '\n';
            std::cout << "Performing sub test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                      << "low:\t" << expect.low << '\t' << result.low << '\n'
                      << "\nhigh:\t" << expect.high << '\t' << result.high << '\n'
                      << "\nexp:\t" << expect.exp << '\t' << result.exp << '\n'
                      << "\nsign:\t" << expect.sign << '\t' << result.sign << '\n';
            return false;
        }
    }
    return true;
}
/*!
 * Test operator*(), return true if passed.
 */
bool numericalMulTest(unsigned int loop) {
    double d;
    for(unsigned int i = 0; i < loop; ++i) {
        d = 1 - 1.4 * d * d;
        double d_a = d * engine();
        MultiScalar a(d_a);

        d = 1 - 1.4 * d * d;
        double d_b = d * engine();
        MultiScalar b(d_b);

        double_extract expect{d_a * d_b};
        double_extract result{double(a * b)};
        if((expect.sign != static_cast<unsigned int>(result.sign))
           || (expect.exp != static_cast<unsigned int>(result.exp))
           || (expect.high != static_cast<unsigned int>(result.high))
           || abs(static_cast<int>(expect.low) - static_cast<int>(result.low)) > 1) {
            std::cout << "Performing mul test " << d_a << " * " << d_b << '\n';
            std::cout << "Performing mul test " "--Failed (" << (i + 1) << '/' << loop << ")\n"
                      << "low:\t" << expect.low << '\t' << result.low << '\n'
                      << "\nhigh:\t" << expect.high << '\t' << result.high << '\n'
                      << "\nexp:\t" << expect.exp << '\t' << result.exp << '\n'
                      << "\nsign:\t" << expect.sign << '\t' << result.sign << '\n';
            return false;
        }
    }
    return true;
}
/*!
 * Test operator/(), return true if passed.
 */
bool numericalDivTest(unsigned int loop) {
    double d;
    for(unsigned int i = 0; i < loop; ++i) {
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
        if((expect.sign != static_cast<unsigned int>(result.sign))
           || (expect.exp != static_cast<unsigned int>(result.exp))
           || (expect.high != static_cast<unsigned int>(result.high))
           || abs(static_cast<int>(expect.low) - static_cast<int>(result.low)) > 1) {
            std::cout << "Performing div test " << d_a << " / " << d_b << '\n';
            std::cout << "Performing div test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                      << "low:\t" << expect.low << '\t' << result.low << '\n'
                      << "high:\t" << expect.high << '\t' << result.high << '\n'
                      << "exp:\t" << expect.exp << '\t' << result.exp << '\n'
                      << "sign:\t" << expect.sign << '\t' << result.sign << '\n';
            return false;
        }
    }
    return true;
}
/*!
 * Tests function toInteger(), return true if passed.
 */
bool toIntegerTest(unsigned int loop) {
    const double max_2 = std::default_random_engine::max() >> 1U;
    for(unsigned int i = 0; i < loop; ++i) {
        //100000 is a arbitrary big number.
        const double temp = (engine() >> 1U) / max_2 * 100000;
        MultiScalar s(temp);
        s.toInteger();
        if(s[0] != static_cast<long>(temp)) { //s and static_cast<long>(temp) must not larger than 2^63.
            std::cout << "toIntegerTest failed: Casting " << temp;
            return false;
        }
    }
    return true;
}

int main(int argc, char** argv) {
    Q_UNUSED(argc)
    Q_UNUSED(argv)
    initPhysica();
    bool result = numericalAddTest(iterateCount)
            && numericalSubTest(iterateCount)
            && numericalMulTest(iterateCount)
            && numericalDivTest(iterateCount)
            && toIntegerTest(iterateCount);
    deInitPhysica();
    return !result;
}