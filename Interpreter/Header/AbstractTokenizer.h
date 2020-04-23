/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_ABSTRACTTOKENIZER_H
#define PHYSICA_ABSTRACTTOKENIZER_H

#include <string>
#include "Token.h"

class AbstractTokenizer {
public:
    AbstractTokenizer() = default;
    AbstractTokenizer(const AbstractTokenizer&) = delete;
    AbstractTokenizer& operator=(const AbstractTokenizer&) = delete;

    const std::string& getBuffer() { return buffer; }
    Token::TokenType getType() { return type; }
protected:
    std::string buffer;
    Token::TokenType type = Token::Null;
    bool stoppable = true;
};

#endif