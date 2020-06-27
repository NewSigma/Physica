/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/Math/Calculus/FunctionTree.h"

namespace Physica::Core {
    /*!
     * Construct a tree which explains the structure of a \class Function.
     * If you want to get a node standing by a constant or a variable, ask it from a \class Function.
     */
    FunctionTree::FunctionTree(MultiScalar (*func)(const MultiScalar&, const MultiScalar&), FunctionTree left, FunctionTree right)
    : func(reinterpret_cast<void*>(func)), left(new FunctionTree(std::move(left))), right(new FunctionTree(std::move(right))) {}

    FunctionTree::FunctionTree(MultiScalar (*func)(const MultiScalar&), FunctionTree f)
    : func(reinterpret_cast<void*>(func)), left(new FunctionTree(std::move(f))), right(nullptr) {}

    FunctionTree::FunctionTree(MultiScalar* constant) //NOLINT No need to initialize left, right.
    : constant(constant), placeHolder(nullptr), func(nullptr) {}

    FunctionTree::FunctionTree(FunctionTree&& func) noexcept : left(func.left), right(func.right), func(func.func) {
        func.left = func.right = nullptr;
    }

    FunctionTree::~FunctionTree() {
        FunctionTree* toDelete = right ? left : nullptr;
        delete right;
        delete toDelete;
    }

    FunctionTree& FunctionTree::operator=(FunctionTree&& f) noexcept {
        left = f.left;
        right = f.right;
        func = f.func;
        f.left = f.right = nullptr;
        return *this;
    }

    MultiScalar FunctionTree::solve() const {
        if(!func)
            return MultiScalar(*constant);
        if(right)
            return reinterpret_cast<MultiScalar (*)(const MultiScalar&, const MultiScalar&)>
            (func)((*left).solve(), (*right).solve());
        return reinterpret_cast<MultiScalar (*)(const MultiScalar&)>(func)((*right).solve());
    }
}