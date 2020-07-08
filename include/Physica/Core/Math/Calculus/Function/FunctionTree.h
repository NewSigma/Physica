/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTIONTREE_H
#define PHYSICA_FUNCTIONTREE_H

#include "Physica/Core/MultiPrecition/Scalar.h"

namespace Physica::Core {
    /*!
     * \class FunctionTree is the intenal data type of a Function.
     */
    class FunctionTree {
    public:
        enum FunctionType {
            Value,
            //functions
            Add,
            Sub,
            Mul,
            Div,
            Square,
            Reciprocal,
            Sqrt,
            Factorial,
            Ln,
            Log,
            Exp,
            Cos,
            Sin,
            Tan,
            Sec,
            Csc,
            Cot,
            ArcCos,
            ArcSin,
            ArcTan,
            ArcSec,
            ArcCsc,
            ArcCot,
            Cosh,
            Sinh,
            Tanh,
            Sech,
            Csch,
            Coth,
            ArcCosh,
            ArcSinh,
            ArcTanh,
            ArcSech,
            ArcCsch,
            ArcCoth,
        };
    private:
        union {
            struct {
                FunctionTree* left;
                FunctionTree* right;
            };
            struct {
                //value must be allocated by \class Function. value must not be deleted by FunctionTree.
                MultiScalar* value{};
                void* placeHolder{};
            };
        };
        //If type equals to nullptr, we use the second struct in the union.
        FunctionType type;
    public:
        FunctionTree(FunctionType type, FunctionTree&& left);
        FunctionTree(FunctionType type, FunctionTree&& left, FunctionTree&& right);
        FunctionTree(FunctionTree&& func) noexcept;
        ~FunctionTree();


        FunctionTree& operator=(FunctionTree&& f) noexcept;
        /* Getters */
        [[nodiscard]] FunctionType getType() const { return type; }
        [[nodiscard]] const FunctionTree* getLeft() const { return type == Value ? nullptr : left; }
        [[nodiscard]] const FunctionTree* getRight() const { return type == Value ? nullptr : right; }
        [[nodiscard]] const MultiScalar* getValue() const { return type == Value ? value : nullptr; }
    private:
        explicit FunctionTree(MultiScalar* value);
        FunctionTree(const FunctionTree& func);

        FunctionTree& operator=(const FunctionTree& func);

        [[nodiscard]] MultiScalar solve() const;
        friend class Function;
        friend class FunctionTreePrinter;
        friend class FunctionDerivative;
    };
}

#endif