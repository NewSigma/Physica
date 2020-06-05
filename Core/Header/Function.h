/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTION_H
#define PHYSICA_FUNCTION_H

#include "Vector.h"
#include "FunctionTree.h"

namespace Physica::Core {
    class Function {
        FunctionTree tree;
        FunctionTree* constantNodes;
        size_t length;
    public:
        Function(FunctionTree tree, size_t variableCount);
        ~Function();

        [[nodiscard]] const FunctionTree& getConstNode(Numerical n, size_t index) const;
    };
}

#endif
