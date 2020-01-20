#include <iostream>
#include <cstring>
#include "../../Header/BasicCalculates.h"
#include "../../Header/Const.h"
#include "../../Header/Solve.h"
//TODO Debug accuracy
/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
extern const Const_1* const_1;
extern const Const_2* const_2;

RealNumberA* reciprocal(const RealNumber* n) {
    return *const_1->ONE / *(RealNumberA*)n;
}

RealNumberA* sqrt(const RealNumber* n) {
    if(!n->sign) {
        std::cout << "[BasicCalculates] Error: Cannot solve the square root of a minus value." << std::endl;
        return nullptr;
    }
    auto MachinePrecision = const_1->MachinePrecision;

    auto copy_n = new RealNumberA;
    *copy_n = *(RealNumberA*)n;
    //Let n < 1 so as to control error.
    char add_power = 0;
    if(copy_n->power > 0) {
        if(copy_n->power % 2 == 0) {
            add_power = char(copy_n->power / 2 + 1);
            copy_n->power = -2;
        }
        else {
            add_power = char((copy_n->power + 1) / 2);
            copy_n->power = -1;
        }
    }
    auto result = getOne();
    RealNumberA* temp;

    //3.33 is the big approximate value of ln(10)/ln(2)
    for(int i = 0; i < 3.33 * MachinePrecision; ++i) {
        temp = *copy_n / *result;
        *result += *temp;
        *result /= *const_1->TWO;
        delete temp;
    }

    result->power += add_power;

    auto error = RealNumberA(new unsigned char[0], 0, -MachinePrecision, true, 1);
    *result += error;

    delete copy_n;
    return result;
}
//TODO not completed: Use gamma function.
RealNumberA* factorial(const RealNumber* n) {
    if(!n->sign) {
        std::cout << "[BasicCalculates] Error: Cannot solve the factorial of a minus value." << std::endl;
        return nullptr;
    }

    RealNumberA* result;
    if(isInteger(n)) {
        result = getOne();
        RealNumber* temp = getOne();

        while(*temp < *n) {
            *temp += *const_1->ONE;
            *result *= *temp;
        }

        delete temp;
        return (RealNumberA*)result;
    }
    return nullptr;
}
//TODO Debug
/*
 * Taylor's formula n.th term: (-1)^n * x^(2n) / (2n!)
 * Here temp_1 = x^(2n), temp_2 = 2n!, rank = 2n
 */
RealNumberA* cos(const RealNumber* n) {
    auto MachinePrecision = const_1->MachinePrecision;
    auto ONE = const_1->ONE;

    auto result = getOne();
    auto square_n = *n * *n;
    auto temp_1 = *n * *n;
    auto temp_2 = getTwo();
    auto rank = getTwo();
    bool sign = false;

    while(true) {
        //Calculate one term of the taylor series.
        auto temp = *temp_1 / *temp_2;
        temp->sign = sign;
        *result += *temp;
        //Here the temp means the criteria of break.
        *temp *= *n;
        *rank += *ONE;
        *temp /= *rank;
        int temp_power = temp->power;
        delete temp;

        if(result->power - temp_power >= MachinePrecision)
            break;
        //Prepare for next calculate.
        sign = !sign;
        *temp_1 *= *square_n;
        *temp_2 *= *rank;
        *rank += *ONE;
        *temp_2 *= *rank;
    }

    auto byte = new unsigned char[MachinePrecision];
    memcpy(byte, result->byte, MachinePrecision * sizeof(char));
    delete[] result->byte;
    result->byte = byte;
    result->length = MachinePrecision;
    result->a = 1;

    delete square_n;
    delete temp_1;
    delete temp_2;
    delete rank;
    return result;
}

RealNumberA* sin(const RealNumber* n) {
    auto MachinePrecision = const_1->MachinePrecision;
    auto ONE = const_1->ONE;

    auto result = getZero();
    auto square_n = *n * *n;
    auto temp_1 = new RealNumberA(n);
    auto temp_2 = getOne();
    auto rank = getOne();
    bool sign = true;

    while(true) {
        //Calculate one term of the taylor series.
        auto temp = *temp_1 / *temp_2;
        temp->sign = sign;
        *result += *temp;
        //Here the temp means the criteria of break.
        *temp *= *n;
        *rank += *ONE;
        *temp /= *rank;
        int temp_power = temp->power;
        delete temp;

        if(result->power - temp_power >= MachinePrecision)
            break;
        //Prepare for next calculate.
        sign = !sign;
        *temp_1 *= *square_n;
        *temp_2 *= *rank;
        *rank += *ONE;
        *temp_2 *= *rank;
    }

    auto byte = new unsigned char[MachinePrecision];
    memcpy(byte, result->byte, MachinePrecision * sizeof(char));
    delete[] result->byte;
    result->byte = byte;
    result->length = MachinePrecision;
    result->a = 1;

    delete square_n;
    delete temp_1;
    delete temp_2;
    delete rank;
    return result;
}

RealNumberA* tan(const RealNumber* n) {
    auto sin_n = sin(n);
    auto cos_n = cos(n);
    auto tan_n = *sin_n / *cos_n;
    delete sin_n;
    delete cos_n;
    return tan_n;
}

RealNumberA* arccos(const RealNumber* n) {
    if(*n == *const_1->ONE)
        return getZero();
    if(*n == *const_1->MINUS_ONE)
        return new RealNumberA(const_1->PI);
    return bisectionMethod(cos, n, const_1->ZERO, const_1->PI, const_1->ONE, const_1->MINUS_ONE);
}

RealNumberA* arcsin(const RealNumber* n) {
    if(*n == *const_1->ONE)
        return new RealNumberA(const_2->PI_DIVIDE_TWO);
    if(*n == *const_1->MINUS_ONE)
        return new RealNumberA(const_2->MINUS_PI_DIVIDE_TWO);
    return bisectionMethod(sin, n, const_2->MINUS_PI_DIVIDE_TWO, const_2->PI_DIVIDE_TWO, const_1->MINUS_ONE, const_1->ONE);
}

RealNumberA* ln(const RealNumber* n) {
    if(!n->sign) {
        std::cout << "[BasicCalculates] Error: Cannot solve the logarithm of a minus value." << std::endl;
        return nullptr;
    }

    auto ONE = const_1->ONE;
    auto MachinePrecision = const_1->MachinePrecision;

    auto result = getZero();
    auto temp_0 = *n + *ONE;
    auto temp_1 = *n - *ONE;
    *temp_1 /= *temp_0;
    auto copy_temp_1 = new RealNumberA(temp_1);
    auto temp_2 = getOne();

    while(true) {
        //Calculate one term of the taylor series.
        auto temp = *temp_1 / *temp_2;
        *result += *temp;
        delete temp;
        //Here the temp means the criteria of break.
        *temp_1 *= *copy_temp_1;
        *temp_2 += *ONE;
        temp = *temp_1 / *temp_2;
        int temp_power = temp->power;
        delete temp;

        if(result->power - temp_power >= MachinePrecision)
            break;
        //Prepare for next calculate.
        *temp_1 *= *copy_temp_1;
        *temp_2 += *ONE;
    }
    *result *= *const_1->TWO;

    auto byte = new unsigned char[MachinePrecision];
    memcpy(byte, result->byte, MachinePrecision * sizeof(char));
    delete[] result->byte;
    result->byte = byte;
    result->length = MachinePrecision;
    result->a = 1;

    delete temp_0;
    delete temp_1;
    delete copy_temp_1;
    delete temp_2;
    return result;
}
//Return log_a n
RealNumberA* log(const RealNumber* n, const RealNumber* a) {
    if(*a == *const_1->ONE) {
        std::cout << "[BasicCalculates] Error: Cannot solve the logarithm of a minus value." << std::endl;
        return nullptr;
    }
    auto ln_n = ln(n);
    auto ln_a = ln(a);
    *ln_n /= *ln_a;
    delete ln_a;
    return ln_n;
}