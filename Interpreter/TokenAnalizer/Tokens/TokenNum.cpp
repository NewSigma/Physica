/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Interpreter/Header/TokenNum.h"
#include "Core/Header/Numerical.h"

using namespace Physica::Interpreter;

TokenNum::TokenNum(const Core::Numerical& num) : Token(TokenType::Numeric), data(new Core::Numerical(num)) {}

TokenNum::~TokenNum() {
    delete data;
}