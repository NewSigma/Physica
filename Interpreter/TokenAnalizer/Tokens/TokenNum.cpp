/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Interpreter/Header/TokenNum.h"
#include "Core/Header/Numerical.h"

TokenNum::TokenNum(const Numerical& num) : Token(TokenType::Numeric), data(new Numerical(num)) {}

TokenNum::~TokenNum() {
    delete data;
}