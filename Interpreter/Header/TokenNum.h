/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOKENNUM_H
#define PHYSICA_TOKENNUM_H

#include "Token.h"

class Numerical;

class TokenNum : public Token {
    Numerical* data;
public:
    explicit TokenNum(const Numerical& num);
    ~TokenNum();
};

#endif