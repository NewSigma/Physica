/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/Function.h"

namespace Physica::Core {
    Function::Function(FunctionTree tree, size_t variableCount)
            : tree(FunctionTree(std::move(tree))), length(variableCount)
            , constantNodes(reinterpret_cast<FunctionTree*>(calloc(variableCount, sizeof(FunctionTree)))) {}

    Function::~Function() {
        auto p = constantNodes;
        for(size_t i = 0; i < length; ++i) {
            if(p)
                p->~FunctionTree();
            ++p;
        }
        free(constantNodes);
    }

    const FunctionTree& Function::getConstNode(Numerical n, size_t index) const {
        const auto p = constantNodes + index;
        if(p == nullptr)
            new (p) FunctionTree(new Numerical(std::move(n)));
        return *p;
    }
}