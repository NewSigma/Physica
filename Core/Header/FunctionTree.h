/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTIONTREE_H
#define PHYSICA_FUNCTIONTREE_H

namespace Physica::Core {
    class Scalar;

    class FunctionTree {
        //A function is a tree.
        union {
            struct {
                FunctionTree* left;
                FunctionTree* right;
            };
            struct {
                Scalar* constant{};
                void* placeHolder{};
            };
        };
        void* func;
    public:
        FunctionTree(Scalar (*func)(const Scalar&, const Scalar&), FunctionTree left, FunctionTree right);
        FunctionTree(Scalar (*func)(const Scalar&), FunctionTree f);
        FunctionTree(const FunctionTree& func) = delete;
        FunctionTree(FunctionTree&& func) noexcept;
        ~FunctionTree();

        FunctionTree& operator=(const FunctionTree& func) = delete;
        FunctionTree& operator=(FunctionTree&& f) noexcept;
    private:
        explicit FunctionTree(Scalar* constant);
        [[nodiscard]] Scalar solve() const;
        friend class Function;
    };
}

#endif