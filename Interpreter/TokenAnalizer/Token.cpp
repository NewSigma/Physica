/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#include <cstring>
#include "Interpreter/Header/Token.h"

Token::Token(const std::string& str, TokenType t) : data(new char[str.length() + 1]), type(t) {
    strcpy(data, str.c_str());
}

Token::~Token() {
    delete[] data;
}