/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTIONTREEPRINTER_H
#define PHYSICA_FUNCTIONTREEPRINTER_H

#include <list>
#include <iosfwd>
#include "Function.h"

namespace Physica::Core {
    class FunctionTreePrinter {
        static const char* rightBranch;
        static const char* leftBranch;
        static const char* straightLine;
        static const char* space;
        static const char* spaces;
        
        std::list<const char*> list;
        const Function& f;
        std::ostream& os;
    public:
        FunctionTreePrinter(const Function& f, std::ostream& os);
        void print() { printImpl(f.getTree(), false); } //true or false leads to the same result.
    private:
        void printImpl(const FunctionTree& functionTree, bool isLeft);
        void printList();
    };
}

#endif