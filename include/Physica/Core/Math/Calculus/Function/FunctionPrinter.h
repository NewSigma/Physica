/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTIONPRINTER_H
#define PHYSICA_FUNCTIONPRINTER_H

#include <iosfwd>
#include "Function.h"

namespace Physica::Core {
    class FunctionPrinter {
        const Function& f;
        std::ostream& os;
    public:
        FunctionPrinter(const Function& f, std::ostream& os);
        void print() { printImpl(f.getTree()); }
    private:
        void printImpl(const FunctionTree& functionTree);
    };
}

#endif