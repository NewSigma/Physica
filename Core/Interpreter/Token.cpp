/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include "../Header/Token.h"

Token::Token(const char* str, TokenType t) : data(str), type(t) {}

Token::~Token() {
    delete[] data;
}