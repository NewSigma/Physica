/*
 * Copyright (c) 2019 NewSigma@163.com. All rights reserved.
 */

#ifndef PHYSICA_TOKENIZERWORD_H
#define PHYSICA_TOKENIZERWORD_H

#include "AbstractTokenizer.h"

class TokenizerWord : public AbstractTokenizer {
    enum BufferState {
        Start,
        Identifier
    } state;
public:
    explicit TokenizerWord(const char*& str);
private:
    bool readChar(char ch);
};

#endif