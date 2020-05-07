/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_TOKENIZER_H
#define PHYSICA_TOKENIZER_H

#include <list>
#include <string>
#include "Token.h"

class Tokenizer {
    enum NumBufferState {
        NumStart,
        Zero,
        Int,
        Float,
        PreExp,
        SignedPreExp,
        Exp
    };

    enum WordBufferState {
        WordStart,
        Identifier
    };
    std::list<Token*> tokens;
    const char* str;
    //Defined to make debug easier.
    int line;
public:
    explicit Tokenizer(const char* str);
    ~Tokenizer();
private:
    void readToken();
    bool readChar(char ch);
    void readNum();
    void readWord();
};

#endif