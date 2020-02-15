#include <iostream>
#include "../../Header/ElementaryFunction.h"
#include "../../Header/Const.h"
#include "../../Header/Solve.h"
//TODO Debug accuracy
/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
extern const Const_1* const_1;
extern const Const_2* const_2;

RealNumber* reciprocal(const RealNumber& n) {
    return *const_1->ONE / n;
}

RealNumber* sqrt_noCheck(const RealNumber& n) {
    if(!n.sign)
        return nullptr;
    auto MachinePrecision = const_1->MachinePrecision;
    auto copy_n = new RealNumber(n);
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

    RealNumber* result = const_1->getOne();
    RealNumber* temp;
    //3.33 is the big approximate value of ln(10)/ln(2)
    for(int i = 0; i < 3.33 * MachinePrecision; ++i) {
        temp = divide(*copy_n, *result);
        *result += *temp;
        *result /= *const_1->TWO;
        delete temp;
    }
    delete copy_n;
    result->power += add_power;
    result->a = 1;

    return result;
}

RealNumber* sqrt(const RealNumber& n) {
    auto result  = sqrt_noCheck(n);
    if(n.a != 0) {
        auto n_error = n.getMinimum();
        auto error = sqrt_noCheck(*n_error);
        if(error == nullptr) {
            delete n_error;
            n_error = n.getMaximum();
            error = sqrt_noCheck(*n_error);
        }
        *error -= *result;
        error->sign = true;

        result->applyError(error);

        delete n_error;
        delete error;
    }
    return result;
}
//TODO not completed: Use gamma function.
RealNumber* factorial(const RealNumber& n) {
    if(!n.sign) {
        std::cout << "[BasicCalculates] Error: Cannot solve the factorial of a minus value." << std::endl;
        return nullptr;
    }

    RealNumber* result;
    if(n.isInteger()) {
        result = const_1->getOne();
        auto temp = const_1->getOne();
        while(*temp < n) {
            *temp += *const_1->ONE;
            *result *= *temp;
        }
        delete temp;
        return result;
    }
    return nullptr;
}
/*
 * Taylor's formula n.th term: (-1)^n * x^(2n) / (2n!)
 * Here temp_1 = x^(2n), temp_2 = 2n!, rank = 2n
 */
RealNumber* cos(const RealNumber& n) {
    auto result = const_1->getOne();
    if(n != *const_1->ZERO) {
        auto MachinePrecision = const_1->MachinePrecision;
        auto ONE = const_1->ONE;

        auto square_n = n * n;
        auto temp_1 = new RealNumber(square_n);
        auto temp_2 = const_1->getTwo();
        auto rank = const_1->getTwo();
        bool sign = false;

        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->sign = sign;
            *result += *temp;
            //Here the temp means the criteria of break.
            *temp *= n;
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
        delete square_n;
        delete temp_1;
        delete temp_2;
        delete rank;
    }
    return result;
}

RealNumber* sin(const RealNumber& n) {
    auto result = const_1->getZero();
    if(n != *const_1->ZERO) {
        auto MachinePrecision = const_1->MachinePrecision;
        auto ONE = const_1->ONE;

        auto square_n = n * n;
        auto temp_1 = new RealNumber(n);
        auto temp_2 = const_1->getOne();
        auto rank = const_1->getOne();
        bool sign = true;

        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->sign = sign;
            *result += *temp;
            //Here the temp means the criteria of break.
            *temp *= n;
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
        delete square_n;
        delete temp_1;
        delete temp_2;
        delete rank;
    }
    return result;
}

RealNumber* tan(const RealNumber& n) {
    auto cos_n = cos(n);
    if(cos_n == nullptr)
        return nullptr;
    auto sin_n = sin(n);
    *sin_n /= *cos_n;
    delete cos_n;
    return sin_n;
}

RealNumber* sec(const RealNumber& n) {
    auto cos_n = cos(n);
    auto result = reciprocal(*cos_n);
    delete cos_n;
    return result;
}

RealNumber* csc(const RealNumber& n) {
    auto sin_n = sin(n);
    auto result = reciprocal(*sin_n);
    delete sin_n;
    return result;
}

RealNumber* cot(const RealNumber& n) {
    auto sin_n = sin(n);
    if(sin_n == nullptr)
        return nullptr;
    auto cos_n = cos(n);
    *cos_n /= *sin_n;
    delete sin_n;
    return cos_n;
}

RealNumber* arccos(const RealNumber& n) {
    return bisectionMethod(cos, n, *const_1->ZERO, *const_2->PI, *const_1->ONE, *const_1->MINUS_ONE);
}

RealNumber* arcsin(const RealNumber& n) {
    return bisectionMethod(sin, n, *const_2->MINUS_PI_DIVIDE_TWO, *const_2->PI_DIVIDE_TWO, *const_1->MINUS_ONE, *const_1->ONE);
}

RealNumber* arctan(const RealNumber& n) {
    auto temp = n * n;
    *temp += *const_1->ONE;
    auto sqrt_temp = sqrt(*temp);
    delete temp;
    
    temp = n / *sqrt_temp;
    delete sqrt_temp;
    
    auto result = arcsin(*temp);
    result->sign = n.sign;
    delete temp;
    return result;
}

RealNumber* arcsec(const RealNumber& n) {
    auto temp = reciprocal(n);
    if(temp == nullptr)
        return temp;
    auto result = arccos(*temp);
    delete temp;
    return result;
}

RealNumber* arccsc(const RealNumber& n) {
    auto temp = reciprocal(n);
    if(temp == nullptr)
        return temp;
    auto result = arcsin(*temp);
    delete temp;
    return result;
}

RealNumber* arccot(const RealNumber& n) {
    auto temp = reciprocal(n);
    if(temp == nullptr)
        return temp;
    auto result = arctan(*temp);
    delete temp;
    return result;
}

RealNumber* ln_noCheck(const RealNumber& n) {
    if(!n.isPositive())
        return nullptr;
    auto result = const_1->getZero();
    if(n != *const_1->ONE) {
        auto temp_0 = n + *const_1->ONE;
        auto temp_1 = n - *const_1->ONE;
        *temp_1 /= *temp_0;
        auto copy_temp_1 = new RealNumber(temp_1);
        auto temp_2 = const_1->getOne();

        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            *result += *temp;
            delete temp;
            //Here the temp means the criteria of break.
            *temp_1 *= *copy_temp_1;
            *temp_2 += *const_1->ONE;
            temp = *temp_1 / *temp_2;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= const_1->MachinePrecision)
                break;
            //Prepare for next calculate.
            *temp_1 *= *copy_temp_1;
            *temp_2 += *const_1->ONE;
        }
        *result *= *const_1->TWO;
        delete temp_0;
        delete temp_1;
        delete copy_temp_1;
        delete temp_2;
    }
    return result;
}

RealNumber* ln(const RealNumber& n) {
    auto result = ln_noCheck(n);
    if(n.a != 0) {
        auto n_error = n.getMinimum();
        auto error = ln_noCheck(*n_error);
        if(error == nullptr) {
            delete n_error;
            n_error = n.getMaximum();
            error = ln_noCheck(*n_error);
        }
        *error -= *result;
        error->sign = true;

        result->applyError(error);

        delete n_error;
        delete error;
    }
    return result;
}
//Return log_a n
RealNumber* log(const RealNumber& n, const RealNumber& a) {
    if(a == *const_1->ONE)
        return nullptr;

    auto ln_n = ln(n);
    if(ln_n == nullptr)
        return nullptr;

    auto ln_a = ln(a);
    if(ln_a == nullptr) {
        delete ln_n;
        return nullptr;
    }
    *ln_n /= *ln_a;
    delete ln_a;
    return ln_n;
}