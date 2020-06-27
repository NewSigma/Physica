/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTIONTREE_H
#define PHYSICA_FUNCTIONTREE_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    class FunctionTree {
        //A function is a tree.
        union {
            struct {
                FunctionTree* left;
                FunctionTree* right;
            };
            struct {
                MultiScalar* constant{};
                void* placeHolder{};
            };
        };
        void* func;
    public:
        FunctionTree(MultiScalar (*func)(const MultiScalar&, const MultiScalar&), FunctionTree left, FunctionTree right);
        FunctionTree(MultiScalar (*func)(const MultiScalar&), FunctionTree f);
        FunctionTree(const FunctionTree& func) = delete;
        FunctionTree(FunctionTree&& func) noexcept;
        ~FunctionTree();

        FunctionTree& operator=(const FunctionTree& func) = delete;
        FunctionTree& operator=(FunctionTree&& f) noexcept;
    private:
        explicit FunctionTree(MultiScalar* constant);
        [[nodiscard]] MultiScalar solve() const;
        friend class Function;
    };
}

#endif