/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOKENWORD_H
#define PHYSICA_TOKENWORD_H

#include "Token.h"

class TokenWord : public Token {
    char* data;
public:
    explicit TokenWord(const char* str);
    TokenWord(const char* str, int len);
    ~TokenWord();
};

#endif