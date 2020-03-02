/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_NUMBER_H
#define PHYSICA_NUMBER_H

//Number is a abstract class.
class Number {
public:
    enum NumberType {
        AbstractNumber,
        RealNumber,
        ComplexNumber,
        RealInf,
        ComplexInf,
        Indeterminate
    };
    NumberType getType();
protected:
    Number() {}
    ~Number() = default;
    Number(const Number& n) = delete;
    NumberType type;
};

#endif
