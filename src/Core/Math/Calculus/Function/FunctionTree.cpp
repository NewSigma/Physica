/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/Math/Calculus/Function/FunctionTree.h"

#include <utility>

namespace Physica::Core {
    /*!
     * Construct a tree which explains the structure of a \class Function.
     * If you want to get a node standing by a value or a variable, ask it from a \class Function.
     */
    FunctionTree::FunctionTree(MultiScalar* value) //NOLINT No need to initialize left, right.
            : value(value), placeHolder(nullptr), type(Value) {}

    FunctionTree::FunctionTree(FunctionType type, FunctionTree&& left)
            : type(type), left(new FunctionTree(std::move(left))), right(nullptr) {}

    FunctionTree::FunctionTree(FunctionType type, FunctionTree&& left, FunctionTree&& right)
    : type(type), left(new FunctionTree(std::move(left))), right(new FunctionTree(std::move(right))) {}
    /*!
     * May be a large cost function. Declared private to avoid incorrect use.
     */
    FunctionTree::FunctionTree(const FunctionTree &func) : type(func.type) {
        //Optimize: may be use operator ? and reinterpret_cast to avoid branches.
        if(func.type == Value) {
            value = func.value;
            placeHolder = nullptr;
        }
        else {
            left = new FunctionTree(*func.left);
            right = new FunctionTree(*func.right);
        }
    }

    FunctionTree::FunctionTree(FunctionTree&& func) noexcept : left(func.left), right(func.right), type(func.type) {
        func.left = func.right = nullptr;
    }

    FunctionTree::~FunctionTree() {
        if(type != Value) {
            delete left;
            delete right;
        }
    }
    /*!
     * May be a large cost function. Declared private to avoid incorrect use.
     */
    FunctionTree& FunctionTree::operator=(const FunctionTree& f) {
        if(this != &f) {
            this->~FunctionTree();
            if(f.type == Value) {
                value = f.value;
                placeHolder = nullptr;
            }
            else {
                left = new FunctionTree(*f.left);
                right = new FunctionTree(*f.right);
            }
            type = f.type;
        }
        return *this;
    }

    FunctionTree& FunctionTree::operator=(FunctionTree&& f) noexcept {
        left = f.left;
        right = f.right;
        type = f.type;
        f.left = f.right = nullptr;
        return *this;
    }

    MultiScalar FunctionTree::solve() const {
        switch(type) {
            case Value:
                return MultiScalar(*value);
            case Add:
                return left->solve() + right->solve();
            case Sub:
                return left->solve() - right->solve();
            case Mul:
                return left->solve() * right->solve();
            case Div:
                return left->solve() / right->solve();
            case Square:
                return square(left->solve());
            case Reciprocal:
                return reciprocal(left->solve());
            case Sqrt:
                return sqrt(left->solve());
            case Factorial:
                return factorial(left->solve());
            case Ln:
                return ln(left->solve());
            case Log:
                return log(left->solve(), right->solve());
            case Exp:
                return exp(left->solve());
            case Cos:
                return cos(left->solve());
            case Sin:
                return sin(left->solve());
            case Tan:
                return tan(left->solve());
            case Sec:
                return sec(left->solve());
            case Csc:
                return csc(left->solve());
            case Cot:
                return cot(left->solve());
            case ArcCos:
                return arccos(left->solve());
            case ArcSin:
                return arcsin(left->solve());
            case ArcTan:
                return arctan(left->solve());
            case ArcSec:
                return arcsec(left->solve());
            case ArcCsc:
                return arccsc(left->solve());
            case ArcCot:
                return arccot(left->solve());
            case Cosh:
                return cosh(left->solve());
            case Sinh:
                return sinh(left->solve());
            case Tanh:
                return tanh(left->solve());
            case Sech:
                return sech(left->solve());
            case Csch:
                return csch(left->solve());
            case Coth:
                return coth(left->solve());
            case ArcCosh:
                return arccosh(left->solve());
            case ArcSinh:
                return arcsinh(left->solve());
            case ArcTanh:
                return arctanh(left->solve());
            case ArcSech:
                return arcsech(left->solve());
            case ArcCsch:
                return arccsch(left->solve());
            case ArcCoth:
                return arccosh(left->solve());
        }
    }
}