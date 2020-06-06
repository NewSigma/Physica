/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "FunctionTree.h"
#include "Core/Header/Numerical.h"

namespace Physica::Core {
    /*!
     * Construct a tree which explains the structure of a \class Function.
     * If you want to get a node standing by a constant or a variable, ask it from a \class Function.
     */
    FunctionTree::FunctionTree(Numerical (*func)(const Numerical&, const Numerical&), FunctionTree left, FunctionTree right)
    : func(reinterpret_cast<void*>(func)), left(new FunctionTree(std::move(left))), right(new FunctionTree(std::move(right))) {}

    FunctionTree::FunctionTree(Numerical (*func)(const Numerical&), FunctionTree f)
    : func(reinterpret_cast<void*>(func)), left(new FunctionTree(std::move(f))), right(nullptr) {}

    FunctionTree::FunctionTree(Numerical* constant) //NOLINT No need to initialize left, right.
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

    Numerical FunctionTree::solve() const {
        if(!func)
            return Numerical(*constant);
        if(right)
            return reinterpret_cast<Numerical (*)(const Numerical&, const Numerical&)>
            (func)((*left).solve(), (*right).solve());
        return reinterpret_cast<Numerical (*)(const Numerical&)>(func)((*right).solve());
    }
}