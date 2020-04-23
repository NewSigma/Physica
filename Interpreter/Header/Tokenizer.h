/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_TOKENIZER_H
#define PHYSICA_TOKENIZER_H

#include <list>
#include <string>
#include "Token.h"

class Tokenizer {
    std::list<Token> tokens;
    const char* str;
    //Defined to make debug easier.
    int line;
public:
    Tokenizer(const char* str);
private:
    void readWord();
    void handleOperator(Token::TokenType single, Token::TokenType pair);
};

#endif