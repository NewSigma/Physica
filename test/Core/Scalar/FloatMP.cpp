/*
 * Copyright 2019-2025 Weibo He.
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
#include <iostream>
#include <random>
#include "Physica/Core/Math/Algebra/LinearAlgebra/Vector/DenseVector.h"
#include "Physica/Core/Scalar/Real.h"
#include "Test.h"

using namespace Physica;
using T = Real<FloatMP>;
constexpr unsigned int iterateCount = 50;
static std::default_random_engine engine(clock());

namespace {
    bool numericalAddTest(unsigned int loop) {
        double d{};
        for (unsigned int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            T a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            T b(d_b);

            double_extract ext{d_a + d_b};
            double_extract result{double(a + b)};
            if ((ext.sign != static_cast<unsigned int>(result.sign))
                || (ext.exp != static_cast<unsigned int>(result.exp))
                || (ext.high != static_cast<unsigned int>(result.high))
                || abs(static_cast<int>(ext.low) - static_cast<int>(result.low)) > 1) {
                std::cout << "Performing add test " << d_a << " + " << d_b << '\n';
                std::cout << "Performing add test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                          << "low:\t" << ext.low << '\t' << result.low << '\n'
                          << "\nhigh:\t" << ext.high << '\t' << result.high << '\n'
                          << "\nexp:\t" << ext.exp << '\t' << result.exp << '\n'
                          << "\nsign:\t" << ext.sign << '\t' << result.sign << '\n';
                return false;
            }
        }
        return true;
    }

    bool numericalSubTest(unsigned int loop) {
        double d{};
        for (unsigned int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            T a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            T b(d_b);

            double_extract ext{d_a - d_b};
            double_extract result{double(a - b)};
            if ((ext.sign != static_cast<unsigned int>(result.sign))
                || (ext.exp != static_cast<unsigned int>(result.exp))
                || (ext.high != static_cast<unsigned int>(result.high))
                || abs(static_cast<int>(ext.low) - static_cast<int>(result.low)) > 1) {
                std::cout << "Performing sub test " << d_a << " - " << d_b << '\n';
                std::cout << "Performing sub test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                          << "low:\t" << ext.low << '\t' << result.low << '\n'
                          << "\nhigh:\t" << ext.high << '\t' << result.high << '\n'
                          << "\nexp:\t" << ext.exp << '\t' << result.exp << '\n'
                          << "\nsign:\t" << ext.sign << '\t' << result.sign << '\n';
                return false;
            }
        }
        return true;
    }

    bool numericalMulTest(unsigned int loop) {
        double d{};
        for (unsigned int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            T a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            T b(d_b);

            double_extract ext{d_a * d_b};
            double_extract result{double(a * b)};
            if ((ext.sign != static_cast<unsigned int>(result.sign))
                || (ext.exp != static_cast<unsigned int>(result.exp))
                || (ext.high != static_cast<unsigned int>(result.high))
                || abs(static_cast<int>(ext.low) - static_cast<int>(result.low)) > 1) {
                std::cout << "Performing mul test " << d_a << " * " << d_b << '\n';
                std::cout << "Performing mul test "
                             "--Failed ("
                          << (i + 1) << '/' << loop << ")\n"
                          << "low:\t" << ext.low << '\t' << result.low << '\n'
                          << "\nhigh:\t" << ext.high << '\t' << result.high << '\n'
                          << "\nexp:\t" << ext.exp << '\t' << result.exp << '\n'
                          << "\nsign:\t" << ext.sign << '\t' << result.sign << '\n';
                return false;
            }
        }
        return true;
    }

    bool numericalDivTest(unsigned int loop) {
        double d{};
        for (unsigned int i = 0; i < loop; ++i) {
            d = 1 - 1.4 * d * d;
            double d_a = d * engine();
            T a(d_a);

            d = 1 - 1.4 * d * d;
            double d_b = d * engine();
            while (d_b == 0) {
                d = 1 - 1.4 * d * d;
                d_b = d * engine();
            }
            T b(d_b);

            double_extract ext{d_a / d_b};
            double_extract result{double(a / b)};
            if ((ext.sign != static_cast<unsigned int>(result.sign))
                || (ext.exp != static_cast<unsigned int>(result.exp))
                || (ext.high != static_cast<unsigned int>(result.high))
                || abs(static_cast<int>(ext.low) - static_cast<int>(result.low)) > 1) {
                std::cout << "Performing div test " << d_a << " / " << d_b << '\n';
                std::cout << "Performing div test " << "--Failed (" << (i + 1) << '/' << loop << ")\n"
                          << "low:\t" << ext.low << '\t' << result.low << '\n'
                          << "high:\t" << ext.high << '\t' << result.high << '\n'
                          << "exp:\t" << ext.exp << '\t' << result.exp << '\n'
                          << "sign:\t" << ext.sign << '\t' << result.sign << '\n';
                return false;
            }
        }
        return true;
    }
}

int main() {
    bool passed = numericalAddTest(iterateCount)
               && numericalSubTest(iterateCount)
               && numericalMulTest(iterateCount)
               && numericalDivTest(iterateCount);
    expect(passed);
    expect(float64(Real<FloatMP>(0.5)) == 0.5);
    {
        // Test that addArrWithArr and addArrWithArrEq round correctlly
        T a({3563280027363695401, 17475862807287258288UL, 9708812670373448218UL, 536}, 4, 0);
        T b({9410232262925970914UL, 16029360673564519792UL, 970881267037344821, 8737931403336103397}, 4, -1);
        expect((a + b)[3] == 537);
    }
    {
        // Test that correctly round if length overflow
        const T a = T({10531239283933347840UL}, -1, -1);
        const T b = T({0, 16275461880519395517UL, 17224891941645282124UL, 11134621108104132433UL}, -4, -1);
        expect((a + b)[3] == 1);
    }
    return 0;
}
