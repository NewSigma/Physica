/*
 * Copyright 2019 WeiBo He.
 *
 * This file is part of Physica.

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
#include <Physica/Physica.h>
#include <Physica/Core/MultiPrecition/Scalar.h>
#include <Physica/Core/MultiPrecition/ScalarImpl/ElementaryFunction.h>

using namespace Physica::Core;

int main(int argc, char** argv) {
    initPhysica();
    std::cout << "Please input a double. Input 0 to stop the test." << '\n';
    double d1, d2;
    std::cin >> d1 >> d2;
    while(d1 != 0) {
        MultiScalar n1(d1);
        MultiScalar n2(d2);
        std::cout << "#######################Basic Functions#######################\n";
        std::cout << "Random: " << randomScalar<MultiPrecision, false>() << '\n';
        std::cout << "Reciprocal: " << reciprocal(n1) << '\n';
        std::cout << "Sqrt: " << sqrt(n1) << '\n';
        if(n1.isInteger())
            std::cout << "Factorial: " << factorial(n1) << '\n';
        std::cout << "ln: " << ln(n1) << '\n';
        std::cout << "log: " << log(n1, n2) << '\n';
        std::cout << "exp: " << exp(n1) << '\n';
        std::cout << "#######################Trigonometric Functions#######################\n";
        std::cout << "cos: " << cos(n1) << '\n';
        std::cout << "sin: " << sin(n1) << '\n';
        std::cout << "tan: " << tan(n1) << '\n';
        std::cout << "sec: " << sec(n1) << '\n';
        std::cout << "csc: " << csc(n1) << '\n';
        std::cout << "cot: " << cot(n1) << '\n';
        std::cout << "arccos: " << arccos(n1) << '\n';
        std::cout << "arcsin: " << arcsin(n1) << '\n';
        std::cout << "arctan: " << arctan(n1) << '\n';
        std::cout << "arcsec: " << arcsec(n1) << '\n';
        std::cout << "arccsc: " << arccsc(n1) << '\n';
        std::cout << "arccot: " << arccot(n1) << '\n';
        std::cout << "#######################Hyperbolic Functions#######################\n";
        std::cout << "cosh: " << cosh(n1) << '\n';
        std::cout << "sinh: " << sinh(n1) << '\n';
        std::cout << "tanh: " << tanh(n1) << '\n';
        std::cout << "sech: " << sech(n1) << '\n';
        std::cout << "csch: " << csch(n1) << '\n';
        std::cout << "coth: " << coth(n1) << '\n';
        std::cout << "arccosh: " << arccosh(n1) << '\n';
        std::cout << "arcsinh: " << arcsinh(n1) << '\n';
        std::cout << "arctanh: " << arctanh(n1) << '\n';
        std::cout << "arcsech: " << arcsech(n1) << '\n';
        std::cout << "arccsch: " << arccsch(n1) << '\n';
        std::cout << "arccoth: " << arccoth(n1) << '\n';

        std::cout << "Waiting for next Input: " << std::endl;
        std::cin >> d1 >> d2;
    }
    deInitPhysica();
    return 0;
}
