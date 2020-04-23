/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOKENIZERNUM_H
#define PHYSICA_TOKENIZERNUM_H

#include "AbstractTokenizer.h"

class TokenizerNum : public AbstractTokenizer {
    enum BufferState {
        Start,
        Zero,
        Int,
        Float,
        PreExp,
        SignedPreExp,
        Exp
    } state;
public:
    explicit TokenizerNum(const char*& str);
private:
    bool readChar(char ch);
};

#endif