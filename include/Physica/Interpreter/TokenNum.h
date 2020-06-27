/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TOKENNUM_H
#define PHYSICA_TOKENNUM_H

#include <Physica/Core/MultiPrecition/Scalar.h>
#include "Token.h"

using Physica::Core::MultiScalar;

namespace Physica::Interpreter {
    class TokenNum : public Token {
        MultiScalar* data;
    public:
        explicit TokenNum(const MultiScalar& num);
        ~TokenNum();
    };
}

#endif