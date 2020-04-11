/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ELEMENTARYFUNCTIONTEST_H
#define PHYSICA_ELEMENTARYFUNCTIONTEST_H

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
        auto random = randomNumerical();
        auto reciprocal_value = reciprocal(n1);
        auto sqrt_value = sqrt(n1);
        Numerical* factorial_value = nullptr;
        if(n1.isInteger())
            factorial_value = factorial(n1);
        auto ln_value = ln(n1);
        auto log_value = log(n1, n2);
        auto exp_value = exp(n1);
        auto pow_value = pow(n1, n2);
        auto cos_value = cos(n1);
        auto sin_value = sin(n1);
        auto tan_value = tan(n1);
        auto sec_value = sec(n1);
        auto csc_value = csc(n1);
        auto cot_value = cot(n1);
        auto arccos_value = arccos(n1);
        auto arcsin_value = arcsin(n1);
        auto arctan_value = arctan(n1);
        auto arcsec_value = arcsec(n1);
        auto arccsc_value = arccsc(n1);
        auto arccot_value = arccot(n1);
        auto cosh_value = cosh(n1);
        auto sinh_value = sinh(n1);
        auto tanh_value = tanh(n1);
        auto sech_value = sech(n1);
        auto csch_value = csch(n1);
        auto coth_value = coth(n1);
        auto arccosh_value = arccosh(n1);
        auto arcsinh_value = arcsinh(n1);
        auto arctanh_value = arctanh(n1);
        auto arcsech_value = arcsech(n1);
        auto arccsch_value = arccsch(n1);
        auto arccoth_value = arccoth(n1);

        std::cout << "#######################Basic Functions#######################\n";
        std::cout << "Random: " << *random << '\n';
        std::cout << "Reciprocal: " << *reciprocal_value << '\n';
        std::cout << "Sqrt: " << *sqrt_value << '\n';
        if(factorial_value != nullptr)
            std::cout << "Factorial: " << *factorial_value << '\n';
        std::cout << "ln: " << *ln_value << '\n';
        std::cout << "log: " << *log_value << '\n';
        std::cout << "exp: " << *exp_value << '\n';
        std::cout << "pow: " << *pow_value << '\n';
        std::cout << "#######################Trigonometric Functions#######################\n";
        std::cout << "cos: " << *cos_value << '\n';
        std::cout << "sin: " << *sin_value << '\n';
        std::cout << "tan: " << *tan_value << '\n';
        std::cout << "sec: " << *sec_value << '\n';
        std::cout << "csc: " << *csc_value << '\n';
        std::cout << "cot: " << *cot_value << '\n';
        std::cout << "arccos: " << *arccos_value << '\n';
        std::cout << "arcsin: " << *arcsin_value << '\n';
        std::cout << "arctan: " << *arctan_value << '\n';
        std::cout << "arcsec: " << *arcsec_value << '\n';
        std::cout << "arccsc: " << *arccsc_value << '\n';
        std::cout << "arccot: " << *arccot_value << '\n';
        std::cout << "#######################Hyperbolic Functions#######################\n";
        std::cout << "cosh: " << *cosh_value << '\n';
        std::cout << "sinh: " << *sinh_value << '\n';
        std::cout << "tanh: " << *tanh_value << '\n';
        std::cout << "sech: " << *sech_value << '\n';
        std::cout << "csch: " << *csch_value << '\n';
        std::cout << "coth: " << *coth_value << '\n';
        std::cout << "arccosh: " << *arccosh_value << '\n';
        std::cout << "arcsinh: " << *arcsinh_value << '\n';
        std::cout << "arctanh: " << *arctanh_value << '\n';
        std::cout << "arcsech: " << *arcsech_value << '\n';
        std::cout << "arccsch: " << *arccsch_value << '\n';
        std::cout << "arccoth: " << *arccoth_value << '\n';

        delete random;
        delete reciprocal_value;
        delete sqrt_value;
        delete factorial_value;
        delete ln_value;
        delete log_value;
        delete exp_value;
        delete pow_value;
        delete cos_value;
        delete sin_value;
        delete tan_value;
        delete sec_value;
        delete csc_value;
        delete cot_value;
        delete arccos_value;
        delete arcsin_value;
        delete arctan_value;
        delete arcsec_value;
        delete arccsc_value;
        delete arccot_value;
        delete cosh_value;
        delete sinh_value;
        delete tanh_value;
        delete sech_value;
        delete csch_value;
        delete coth_value;
        delete arccosh_value;
        delete arcsinh_value;
        delete arctanh_value;
        delete arcsech_value;
        delete arccsch_value;
        delete arccoth_value;
        std::cout << "Waiting for next Input: " << std::endl;
        std::cin >> d1 >> d2;
    }
}

#endif //PHYSICA_ELEMENTARYFUNCTIONTEST_H