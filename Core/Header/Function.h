/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTION_H
#define PHYSICA_FUNCTION_H

#include <utility>

#include "Vector.h"
#include "FunctionTree.h"

namespace Physica::Core {
    class Function {
        FunctionTree tree;
        FunctionTree* constantNodes;
        size_t length;
    public:
        Function(FunctionTree tree, size_t variableCount);
        Function(const Function& func) = delete;
        Function(Function&& func) noexcept;
        ~Function();

        Function& operator=(const Function& func) = delete;
        Function& operator=(Function&& func) noexcept;
        [[nodiscard]] inline Numerical& operator[](size_t index);
        [[nodiscard]] inline const Numerical& operator[](size_t index) const;
        [[nodiscard]] Numerical operator()(Numerical n);
        [[nodiscard]] Numerical operator()(Numerical n1, Numerical n2);
        [[nodiscard]] Numerical operator()(Numerical n1, Numerical n2, Numerical n3);

        [[nodiscard]] const FunctionTree& getConstNode(size_t index) const { return constantNodes[index]; }
        [[nodiscard]] Numerical solve() const { return tree.solve(); }
    };
    /* Inline implementations */
    inline Numerical& Function::operator[](size_t index) {
        Q_ASSERT(index < length);
        return *constantNodes[index].constant;
    }

    inline const Numerical& Function::operator[](size_t index) const {
        Q_ASSERT(index < length);
        return *constantNodes[index].constant;
    }
}

#endif
