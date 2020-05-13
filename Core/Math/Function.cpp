/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Function.h"

namespace Physica::Core {
    Function::Function(Numerical (*func)(const Numerical&, const Numerical&), Function* left, Function* right)
    : func(reinterpret_cast<void*>(func)), left(left), right(right), c(nullptr) {}

    Function::Function(const Numerical& n) : func(nullptr), left(nullptr), right(nullptr), c(new Numerical(n)) {}

    Function::~Function() {
        delete left;
        delete right;
        delete c;
    }
    /*
     * Warning:
     * Here we may cast (Numerical (*)(const Numerical&)) to (Numerical (*)(const Numerical&, const Numerical&)),
     * this operation will use one more register and may not cause other effects (On G++ 7.5.0)
     * and will execute *nullptr while we will not use the value.
     *
     * This function is a test edition for performance, if we have better implementation we may remove it.
     */
    Numerical Function::operator()(const Numerical& n) {
        if(left) {
            if(right)
                return reinterpret_cast<Numerical (*)(const Numerical&, const Numerical&)>(func)((*left)(n), (*right)(n));
            else
                return reinterpret_cast<Numerical (*)(const Numerical&, const Numerical&)>(func)((*left)(n), *c);
        }
        else {
            if(right) {
                if(c)
                    return reinterpret_cast<Numerical (*)(const Numerical&, const Numerical&)>(func)(*c, (*right)(n));
                return reinterpret_cast<Numerical (*)(const Numerical&)>(func)((*right)(n));
            }
        }
        return Numerical(c);
    }
}