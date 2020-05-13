/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTION_H
#define PHYSICA_FUNCTION_H

#include "Numerical.h"

namespace Physica::Core {
    class Function {
        //A function is a tree.
        Function* left;
        Function* right;
        void* func;
        Numerical* c;
    public:
        Function(Numerical (*func)(const Numerical&, const Numerical&), Function* left, Function* right);
        explicit Function(const Numerical& n);
        ~Function();
        Numerical operator()(const Numerical& n);
    };
}

#endif