/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOKENNUM_H
#define PHYSICA_TOKENNUM_H

#include "Token.h"

namespace Physica::Core {
    class Numerical;
}

using Physica::Core::Numerical;

namespace Physica::Interpreter {
    class TokenNum : public Token {
        Core::Numerical* data;
    public:
        explicit TokenNum(const Numerical& num);
        ~TokenNum();
    };
}

#endif