/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Core/Header/Function.h"

namespace Physica::Core {
    Function::Function(FunctionTree tree, size_t variableCount)
            : tree(FunctionTree(std::move(tree))), length(variableCount)
            , constantNodes(reinterpret_cast<FunctionTree*>(calloc(variableCount, sizeof(FunctionTree)))) {
        for(size_t i = 0; i < variableCount; ++i)
            new (constantNodes + i) FunctionTree(new Numerical());
    }

    Function::Function(Function&& func) noexcept
            : tree(std::move(func.tree)), constantNodes(func.constantNodes), length(func.length) {
        func.constantNodes = nullptr;
        func.length = 0;
    }

    Function::~Function() {
        for(size_t i = 0; i < length; ++i)
            delete (constantNodes + i)->constant;
        free(constantNodes);
    }

    Function& Function::operator=(Function&& func) noexcept {
        tree = std::move(func.tree);
        constantNodes = func.constantNodes;
        func.constantNodes = nullptr;
        length = func.length;
        func.length = 0;
        return *this;
    }

    Numerical Function::operator()(Numerical n) {
        (*this)[0] = std::move(n);
        return solve();
    }

    Numerical Function::operator()(Numerical n1, Numerical n2) {
        (*this)[0] = std::move(n1);
        (*this)[1] = std::move(n2);
        return solve();
    }

    Numerical Function::operator()(Numerical n1, Numerical n2, Numerical n3) {
        (*this)[0] = std::move(n1);
        (*this)[1] = std::move(n2);
        (*this)[2] = std::move(n3);
        return solve();
    }
}