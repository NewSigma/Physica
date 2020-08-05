/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TREEFUNCTIONPRINTER_H
#define PHYSICA_TREEFUNCTIONPRINTER_H

#include <list>
#include <iosfwd>
#include "TreeFunctionData.h"

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    class TreeFunction;

    namespace {
        class TreeFunctionPrinterPrivate {
        protected:
            static const char* rightBranch;
            static const char* leftBranch;
            static const char* straightLine;
            static const char* space;
            static const char* spaces;
        };
    }

    template<ScalarType type, bool errorTrack>
    class TreeFunctionPrinter : private TreeFunctionPrinterPrivate {
        std::list<const char*> list;
        const TreeFunction<type, errorTrack>& f;
        std::ostream& os;
    public:
        TreeFunctionPrinter(const TreeFunction<type, errorTrack>& f, std::ostream& os);
        void print() { printImpl(f.getTree(), false); } //true or false leads to the same result.
    private:
        void printImpl(const TreeFunctionData<type, errorTrack>& functionTree, bool isLeft);
        void printList();
    };
}

#include "Physica/Core/Math/Calculus/Function/TreeFunction/TreeFunctionImpl/TreeFunctionDataPrinterImpl.h"

#endif