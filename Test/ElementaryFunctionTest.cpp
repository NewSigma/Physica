/*
  Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <iostream>
#include <Numerical.h>
#include <ElementaryFunction.h>

void elementary_function_test() {
    std::cout << "Please input a double. Input 0 to stop the test." << '\n';
    double d1 = 1, d2;
    std::cin >> d1 >> d2;
    while(d1 != 0) {
        Numerical n1(d1);
        Numerical n2(d2);
        Numerical random = randomNumerical();
        Numerical reciprocal_value = reciprocal(n1);
        Numerical sqrt_value = sqrt(n1);
        Numerical factorial_value = getZero();
        if(n1.isInteger())
            factorial_value = factorial(n1);
        Numerical ln_value = ln(n1);
        Numerical log_value = log(n1, n2);
        Numerical exp_value = exp(n1);
        Numerical pow_value = pow(n1, n2);
        Numerical cos_value = cos(n1);
        Numerical sin_value = sin(n1);
        Numerical tan_value = tan(n1);
        Numerical sec_value = sec(n1);
        Numerical csc_value = csc(n1);
        Numerical cot_value = cot(n1);
        Numerical arccos_value = arccos(n1);
        Numerical arcsin_value = arcsin(n1);
        Numerical arctan_value = arctan(n1);
        Numerical arcsec_value = arcsec(n1);
        Numerical arccsc_value = arccsc(n1);
        Numerical arccot_value = arccot(n1);
        Numerical cosh_value = cosh(n1);
        Numerical sinh_value = sinh(n1);
        Numerical tanh_value = tanh(n1);
        Numerical sech_value = sech(n1);
        Numerical csch_value = csch(n1);
        Numerical coth_value = coth(n1);
        Numerical arccosh_value = arccosh(n1);
        Numerical arcsinh_value = arcsinh(n1);
        Numerical arctanh_value = arctanh(n1);
        Numerical arcsech_value = arcsech(n1);
        Numerical arccsch_value = arccsch(n1);
        Numerical arccoth_value = arccoth(n1);

        std::cout << "#######################Basic Functions#######################\n";
        std::cout << "Random: " << random << '\n';
        std::cout << "Reciprocal: " << reciprocal_value << '\n';
        std::cout << "Sqrt: " << sqrt_value << '\n';
        std::cout << "Factorial: " << factorial_value << '\n';
        std::cout << "ln: " << ln_value << '\n';
        std::cout << "log: " << log_value << '\n';
        std::cout << "exp: " << exp_value << '\n';
        std::cout << "pow: " << pow_value << '\n';
        std::cout << "#######################Trigonometric Functions#######################\n";
        std::cout << "cos: " << cos_value << '\n';
        std::cout << "sin: " << sin_value << '\n';
        std::cout << "tan: " << tan_value << '\n';
        std::cout << "sec: " << sec_value << '\n';
        std::cout << "csc: " << csc_value << '\n';
        std::cout << "cot: " << cot_value << '\n';
        std::cout << "arccos: " << arccos_value << '\n';
        std::cout << "arcsin: " << arcsin_value << '\n';
        std::cout << "arctan: " << arctan_value << '\n';
        std::cout << "arcsec: " << arcsec_value << '\n';
        std::cout << "arccsc: " << arccsc_value << '\n';
        std::cout << "arccot: " << arccot_value << '\n';
        std::cout << "#######################Hyperbolic Functions#######################\n";
        std::cout << "cosh: " << cosh_value << '\n';
        std::cout << "sinh: " << sinh_value << '\n';
        std::cout << "tanh: " << tanh_value << '\n';
        std::cout << "sech: " << sech_value << '\n';
        std::cout << "csch: " << csch_value << '\n';
        std::cout << "coth: " << coth_value << '\n';
        std::cout << "arccosh: " << arccosh_value << '\n';
        std::cout << "arcsinh: " << arcsinh_value << '\n';
        std::cout << "arctanh: " << arctanh_value << '\n';
        std::cout << "arcsech: " << arcsech_value << '\n';
        std::cout << "arccsch: " << arccsch_value << '\n';
        std::cout << "arccoth: " << arccoth_value << '\n';

        std::cout << "Waiting for next Input: " << std::endl;
        std::cin >> d1 >> d2;
    }
}