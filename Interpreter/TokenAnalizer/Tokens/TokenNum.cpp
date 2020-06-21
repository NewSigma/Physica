/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Interpreter/Header/TokenNum.h"
#include "Core/Header/Scalar.h"

using namespace Physica::Interpreter;

TokenNum::TokenNum(const Core::Scalar& num) : Token(TokenType::Numeric), data(new Core::Scalar(num)) {}

TokenNum::~TokenNum() {
    delete data;
}