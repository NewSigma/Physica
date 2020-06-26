/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/Math/Calculus/FunctionTree.h"
#include "Physica/Core/Scalar.h"

namespace Physica::Core {
    /*!
     * Construct a tree which explains the structure of a \class Function.
     * If you want to get a node standing by a constant or a variable, ask it from a \class Function.
     */
    FunctionTree::FunctionTree(Scalar (*func)(const Scalar&, const Scalar&), FunctionTree left, FunctionTree right)
    : func(reinterpret_cast<void*>(func)), left(new FunctionTree(std::move(left))), right(new FunctionTree(std::move(right))) {}

    FunctionTree::FunctionTree(Scalar (*func)(const Scalar&), FunctionTree f)
    : func(reinterpret_cast<void*>(func)), left(new FunctionTree(std::move(f))), right(nullptr) {}

    FunctionTree::FunctionTree(Scalar* constant) //NOLINT No need to initialize left, right.
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

    Scalar FunctionTree::solve() const {
        if(!func)
            return Scalar(*constant);
        if(right)
            return reinterpret_cast<Scalar (*)(const Scalar&, const Scalar&)>
            (func)((*left).solve(), (*right).solve());
        return reinterpret_cast<Scalar (*)(const Scalar&)>(func)((*right).solve());
    }
}