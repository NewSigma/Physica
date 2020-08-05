/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTIONPRINTER_H
#define PHYSICA_FUNCTIONPRINTER_H

#include <iosfwd>
#include "TreeFunction.h"

namespace Physica::Core {
    template<ScalarType type = MultiPrecision, bool errorTrack = true>
    class FunctionPrinter {
        const TreeFunction<type, errorTrack>& f;
        std::ostream& os;
    public:
        FunctionPrinter(const TreeFunction<type, errorTrack>& f, std::ostream& os);
        void print() { printImpl(f.getTree()); }
    private:
        void printImpl(const TreeFunctionData<type, errorTrack>& functionTree);
    };
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/FunctionPrinterImpl.h"

#endif