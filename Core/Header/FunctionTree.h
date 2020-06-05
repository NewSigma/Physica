/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTIONTREE_H
#define PHYSICA_FUNCTIONTREE_H

namespace Physica::Core {
    class Numerical;

    class FunctionTree {
        //A function is a tree.
        union {
            struct {
                FunctionTree* left;
                FunctionTree* right;
            };
            struct {
                const Numerical* constant{};
                void* placeHolder{};
            };
        };
        void* func;
    public:
        FunctionTree(Numerical (*func)(const Numerical&, const Numerical&), FunctionTree left, FunctionTree right);
        FunctionTree(Numerical (*func)(const Numerical&), FunctionTree f);
        FunctionTree(const FunctionTree& func) = delete;
        FunctionTree(FunctionTree&& func) noexcept;
        ~FunctionTree();

        FunctionTree& operator=(const FunctionTree& func) = delete;
        FunctionTree& operator=(FunctionTree&& f) noexcept;

        [[nodiscard]] Numerical solve() const;
    private:
        explicit FunctionTree(const Numerical* constant);
        friend class Function;
    };
}

#endif