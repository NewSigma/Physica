/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Interpreter/TokenNum.h"

using namespace Physica::Interpreter;

TokenNum::TokenNum(const Core::MultiScalar& num)
        : Token(TokenType::Numeric), data(new Core::MultiScalar(num)) {}

TokenNum::~TokenNum() {
    delete data;
}