/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <QtCore/qlogging.h>
#include "Solve.h"
#include "ElementaryFunction.h"
#include "CalcBasic.h"
#include "Numerical.h"

extern const MathConst* mathConst;
extern const unsigned char UnitByte;

//Return a real number between 0 and 1.
Numerical* randomNumerical() {
    srand(clock());
    srand(random());
    auto result = new Numerical((double)random());
    *result /= basicConst->getR_MAX();

    return result;
}
//Return a real number lowerBound and upperBound.
Numerical* randomNumerical(Numerical* lowerBound, Numerical* upperBound) {
    Numerical* random = randomNumerical();
    auto result = *lowerBound - *upperBound;
    *random *= *random;
    *result += *lowerBound;
    delete random;

    return result;
}

Numerical* reciprocal(const Numerical& n) {
    return basicConst->get_1() / n;
}
/*
 * *_light functions do not consider the error caused by a. For example, sqrt_light does not calculate
 * (sqrt(n + a) - sqrt(n)) for error.
 */
Numerical* sqrt_light(const Numerical& n) {
    if(n.length < 0)
        return nullptr;
    auto GlobalPrecision = basicConst->getGlobalPrecision();
    auto copy_n = new Numerical(n);
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

    Numerical* result = getOne();
    Numerical* temp;
    //3.33 is the big approximate value of ln(10)/ln(2)
    for(int i = 0; i < CHAR_BIT * UnitByte * GlobalPrecision; ++i) {
        temp = divide(*copy_n, *result);
        *result += *temp;
        *result /= basicConst->get_2();
        delete temp;
    }
    delete copy_n;
    result->power += add_power;
    result->a = 1;

    return result;
}

Numerical* sqrt(const Numerical& n) {
    auto result  = sqrt_light(n);
    if(n.a != 0) {
        auto n_error = n.getMinimum();
        auto error = sqrt_light(*n_error);
        if(error == nullptr) {
            delete n_error;
            n_error = n.getMaximum();
            error = sqrt_light(*n_error);
        }
        *error -= *result;
        error->length = (signed char)error->getSize();

        result->applyError(error);

        delete n_error;
        delete error;
    }
    return result;
}
//TODO not completed: Use gamma function.
Numerical* factorial(const Numerical& n) {
    if(n.length < 0) {
        qCritical("Cannot solve the factorial of a minus value.");
        return nullptr;
    }

    Numerical* result;
    if(n.isInteger()) {
        result = getOne();
        auto temp = getOne();
        while(*temp < n) {
            *temp += basicConst->get_1();
            *result *= *temp;
        }
        delete temp;
        return result;
    }
    return nullptr;
}

Numerical* ln_light(const Numerical& n) {
    if(!n.isPositive())
        return nullptr;
    auto result = getZero();
    if(n != basicConst->get_1()) {
        auto temp_0 = add(n, basicConst->get_1());
        auto temp_1 = subtract(n, basicConst->get_1());
        *temp_1 /= *temp_0;
        auto copy_temp_1 = new Numerical(temp_1);
        auto temp_2 = getOne();

        copy_temp_1->a = temp_1->a = 0;
        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->a = 0;
            *result += *temp;
            delete temp;
            //Here the temp means the criteria of break.
            *temp_1 *= *copy_temp_1;
            *temp_2 += basicConst->get_1();
            temp = *temp_1 / *temp_2;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= basicConst->getGlobalPrecision())
                break;
            //Prepare for next calculate.
            *temp_1 *= *copy_temp_1;
            *temp_2 += basicConst->get_1();
        }
        *result *= basicConst->get_2();
        delete temp_0;
        delete temp_1;
        delete copy_temp_1;
        delete temp_2;
    }
    return result;
}

Numerical* ln(const Numerical& n) {
    auto result = ln_light(n);
    if(n.a != 0) {
        auto n_error = n.getMinimum();
        auto error = ln_light(*n_error);
        if(error == nullptr) {
            delete n_error;
            n_error = n.getMaximum();
            error = ln_light(*n_error);
        }
        *error -= *result;
        error->length = (signed char)error->getSize();

        result->applyError(error);

        delete n_error;
        delete error;
    }
    return result;
}
//Return log_a n
Numerical* log(const Numerical& n, const Numerical& a) {
    if(a == basicConst->get_1())
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

Numerical* exp(const Numerical& n) {
    auto result = getOne();
    auto temp = new Numerical(n);
    auto rank = getOne();
    while(true) {
        *temp /= *rank;
        if(*temp < basicConst->getExpectedRelativeError())
            break;
        *result += *temp;
        *temp *= n;
        *rank += basicConst->get_1();
    }
    return result;
}

Numerical* pow(const Numerical& n, const Numerical& a) {
    auto result = ln(n);
    *result *= a;
    *result << *exp(*result);
    return result;
}
/*
 * Taylor's formula n.th term: (-1)^n * x^(2n) / (2n!)
 * Here temp_1 = x^(2n), temp_2 = 2n!, rank = 2n
 */
Numerical* cos(const Numerical& n) {
    auto result = getOne();
    if(n != basicConst->get_0()) {
        auto& ONE = basicConst->get_1();

        auto square_n = n * n;
        auto temp_1 = new Numerical(square_n);
        auto temp_2 = getTwo();
        auto rank = getTwo();

        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->length = (signed char)-temp->length;
            *result += *temp;
            //Here the temp means the criteria of break.
            *temp *= n;
            *rank += ONE;
            *temp /= *rank;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= basicConst->getGlobalPrecision())
                break;
            //Prepare for next calculate.
            *temp_1 *= *square_n;
            *temp_2 *= *rank;
            *rank += ONE;
            *temp_2 *= *rank;
        }
        delete square_n;
        delete temp_1;
        delete temp_2;
        delete rank;
    }
    return result;
}

Numerical* sin(const Numerical& n) {
    auto result = getZero();
    if(n != basicConst->get_0()) {
        auto& ONE = basicConst->get_1();

        auto square_n = n * n;
        auto temp_1 = new Numerical(n);
        auto temp_2 = getOne();
        auto rank = getOne();

        while(true) {
            //Calculate one term of the taylor series.
            auto temp = *temp_1 / *temp_2;
            temp->length = (signed char)-temp->length;
            *result += *temp;
            //Here the temp means the criteria of break.
            *temp *= n;
            *rank += ONE;
            *temp /= *rank;
            int temp_power = temp->power;
            delete temp;

            if(result->power - temp_power >= basicConst->getGlobalPrecision())
                break;
            //Prepare for next calculate.
            *temp_1 *= *square_n;
            *temp_2 *= *rank;
            *rank += ONE;
            *temp_2 *= *rank;
        }
        delete square_n;
        delete temp_1;
        delete temp_2;
        delete rank;
    }
    return result;
}

Numerical* tan(const Numerical& n) {
    auto cos_n = cos(n);
    if(cos_n == nullptr)
        return nullptr;
    auto sin_n = sin(n);
    *sin_n /= *cos_n;
    delete cos_n;
    return sin_n;
}

Numerical* sec(const Numerical& n) {
    auto cos_n = cos(n);
    auto result = reciprocal(*cos_n);
    delete cos_n;
    return result;
}

Numerical* csc(const Numerical& n) {
    auto sin_n = sin(n);
    auto result = reciprocal(*sin_n);
    delete sin_n;
    return result;
}

Numerical* cot(const Numerical& n) {
    auto sin_n = sin(n);
    if(sin_n == nullptr)
        return nullptr;
    auto cos_n = cos(n);
    *cos_n /= *sin_n;
    delete sin_n;
    return cos_n;
}

Numerical* arccos(const Numerical& n) {
    return bisectionMethod(cos, n, basicConst->get_0(), mathConst->getPI(), basicConst->get_1(), basicConst->getMinus_1());
}

Numerical* arcsin(const Numerical& n) {
    return bisectionMethod(sin, n, mathConst->getMinus_PI_2(), mathConst->getPI_2(), basicConst->getMinus_1(), basicConst->get_1());
}

Numerical* arctan(const Numerical& n) {
    auto temp = n * n;
    *temp += basicConst->get_1();
    auto sqrt_temp = sqrt(*temp);
    delete temp;

    temp = n / *sqrt_temp;
    delete sqrt_temp;

    auto result = arcsin(*temp);
    if((result->length ^ n.length) < 0) // NOLINT(hicpp-signed-bitwise)
        result->length = (signed char)-result->length;

    delete temp;
    return result;
}

Numerical* arcsec(const Numerical& n) {
    auto temp = reciprocal(n);
    if(temp == nullptr)
        return temp;
    auto result = arccos(*temp);
    delete temp;
    return result;
}

Numerical* arccsc(const Numerical& n) {
    auto temp = reciprocal(n);
    if(temp == nullptr)
        return temp;
    auto result = arcsin(*temp);
    delete temp;
    return result;
}

Numerical* arccot(const Numerical& n) {
    auto temp = reciprocal(n);
    if(temp == nullptr)
        return temp;
    auto result = arctan(*temp);
    delete temp;
    return result;
}

Numerical* cosh(const Numerical& n) {
    auto result = exp(n);
    auto temp = reciprocal(*result);
    *result += *temp;
    *result /= basicConst->get_2();
    delete temp;
    return result;
}

Numerical* sinh(const Numerical& n) {
    auto result = exp(n);
    auto temp = reciprocal(*result);
    *result -= *temp;
    *result /= basicConst->get_2();
    delete temp;
    return result;
}

Numerical* tanh(const Numerical& n) {
    auto result = exp(n);
    auto temp = reciprocal(*result);
    auto temp1 = *result + *temp;
    *result -= *temp;
    *result /= *temp1;
    delete temp;
    delete temp1;
    return result;
}

Numerical* sech(const Numerical& n) {
    auto result = getTwo();
    auto temp = exp(n);
    auto temp1 = reciprocal(*temp);
    *temp += *temp1;
    *result /= *temp;
    delete temp;
    delete temp1;
    return result;
}

Numerical* csch(const Numerical& n) {
    auto result = getTwo();
    auto temp = exp(n);
    auto temp1 = reciprocal(*temp);
    *temp -= *temp1;
    *result /= *temp;
    delete temp;
    delete temp1;
    return result;
}

Numerical* coth(const Numerical& n) {
    auto result = exp(n);
    auto temp = reciprocal(*result);
    auto temp1 = *result - *temp;
    *result += *temp;
    *result /= *temp1;
    delete temp;
    delete temp1;
    return result;
}

Numerical* arccosh(const Numerical& n) {
    auto temp = n * n;
    *temp -= basicConst->get_1();
    *temp << *sqrt(*temp);
    *temp += n;
    auto result = ln(*temp);
    delete temp;
    return result;
}

Numerical* arcsinh(const Numerical& n) {
    auto temp = n * n;
    *temp += basicConst->get_1();
    *temp << *sqrt(*temp);
    *temp += n;
    auto result = ln(*temp);
    delete temp;
    return result;
}

Numerical* arctanh(const Numerical& n) {
    auto result = basicConst->get_1() + n;
    auto temp = basicConst->get_1() - n;
    *result /= *temp;
    *result << *ln(*result);
    *result /= basicConst->get_2();
    delete temp;
    return result;
}

Numerical* arcsech(const Numerical& n) {
    auto temp = reciprocal(n);
    auto result = arccosh(*temp);
    delete temp;
    return result;
}

Numerical* arccsch(const Numerical& n) {
    auto temp = reciprocal(n);
    auto result = arcsinh(*temp);
    delete temp;
    return result;
}

Numerical* arccoth(const Numerical& n) {
    auto result = n + basicConst->get_1();
    auto temp = n - basicConst->get_1();
    *result /= *temp;
    *result << *ln(*result);
    *result /= basicConst->get_2();
    delete temp;
    return result;
}