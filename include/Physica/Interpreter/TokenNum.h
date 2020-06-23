/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOKENNUM_H
#define PHYSICA_TOKENNUM_H

#include "Token.h"

namespace Physica::Core {
    class Scalar;
}

using Physica::Core::Scalar;

namespace Physica::Interpreter {
    class TokenNum : public Token {
        Core::Scalar* data;
    public:
        explicit TokenNum(const Scalar& num);
        ~TokenNum();
    };
}

#endif