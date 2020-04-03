/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_TOKENIZER_H
#define PHYSICA_TOKENIZER_H

#include <list>
#include <string>
#include "Token.h"

class Tokenizer {
    enum BufferState {
        Start,
        PreInt,
        Int,
        Float
    } state;
    std::list<Token> tokens;
    std::string buffer;
public:
    Tokenizer(const char* str);
    ~Tokenizer() = default;

    inline bool isNonZeroNum(const char& ch) { return '1' <= ch && ch <= '9'; }
private:
    void readChar(const char& ch);
    void resetBuffer(Token::TokenType type);
    Token::TokenType getCurrentType();
};

#endif