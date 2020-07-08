/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/Math/Calculus/Function/FunctionDerivative.h"
#include "Physica/Core/Math/Calculus/Function/Function.h"

namespace Physica::Core {
    FunctionDerivative::FunctionDerivative(const Function& f, size_t index) : f(f), index(index) {}
    /*!
     * Return the partial derivative of the given function.
     */
    Function FunctionDerivative::derivative() const {
        const auto variablesLength = f.getVariablesLength();
        const auto constantsLength = f.getConstantsLength() + 4;
        Function result(derivativeTree(f.getTree()), variablesLength, constantsLength);
        for(size_t i = 0; i < variablesLength; ++i)
            result.setVariable(f.getVariable(i), i);
        //Generate some constant values.
        result.setConstant(MultiScalar::getZero(), 0);
        result.setConstant(MultiScalar::getOne(), 1);
        result.setConstant(MultiScalar::getTwo(), 2);
        result.setConstant(-MultiScalar::getOne(), 3);
        for(size_t i = 3; i < constantsLength; ++i)
            result.setConstant(f.getConstant(i), i);
        return result;
    }
    /*!
     * Return the partial derivative of the given tree.
     *
     * Issue:
     * The time complexity may be very large.
     *
     * Optimize:
     * Derivative of value may return 0 or 1. Make use of them.
     */
    FunctionTree FunctionDerivative::derivativeTree(const FunctionTree& tree) const {
        const FunctionTree::FunctionType type = tree.getType();
        switch(type) {
            case FunctionTree::Value:
                return index == f.getVariablePos(tree) ? f.getConstantNode(1) : f.getConstantNode(0);
            case FunctionTree::Add:
            case FunctionTree::Sub:
                return FunctionTree(type
                        , derivativeTree(*tree.getLeft()), derivativeTree(*tree.getLeft()));
            case FunctionTree::Mul:
                return FunctionTree(FunctionTree::Add
                        , FunctionTree(type, derivativeTree(*tree.getLeft()), FunctionTree(*tree.getRight()))
                        , FunctionTree(type, derivativeTree(*tree.getRight()), FunctionTree(*tree.getLeft())));
            case FunctionTree::Div: {
                FunctionTree f1(FunctionTree::Mul, FunctionTree(*tree.getRight()), derivativeTree(*tree.getLeft()));
                FunctionTree f2(FunctionTree::Mul, FunctionTree(*tree.getLeft()), derivativeTree(*tree.getRight()));
                FunctionTree f1_f2(FunctionTree::Sub, std::move(f1), std::move(f2));
                FunctionTree f2_2(FunctionTree::Square, FunctionTree(*tree.getRight()));
                return FunctionTree(FunctionTree::Div, std::move(f1_f2), std::move(f2_2));
            }
            case FunctionTree::Square: {
                return FunctionTree(FunctionTree::Mul, derivativeTree(*tree.getLeft())
                        , FunctionTree(FunctionTree::Mul, f.getConstantNode(2), FunctionTree(*tree.getLeft())));
            }
            case FunctionTree::Reciprocal: {
                FunctionTree s(FunctionTree::Square, FunctionTree(*tree.getLeft()));
                return FunctionTree(FunctionTree::Mul
                        , f.getConstantNode(3)
                        , FunctionTree(FunctionTree::Div, derivativeTree(*tree.getLeft()), std::move(s)));
            }
            case FunctionTree::Sqrt: {
                FunctionTree s(FunctionTree::Sqrt, FunctionTree(*tree.getLeft()));
                FunctionTree r(FunctionTree::Div, derivativeTree(*tree.getLeft()), std::move(s));
                return FunctionTree(FunctionTree::Div, std::move(r), f.getConstantNode(2));
            }
            case FunctionTree::Factorial:
                //Error
                return f.getConstantNode(0);
            case FunctionTree::Ln:
                return FunctionTree(FunctionTree::Div
                        , derivativeTree(*tree.getLeft())
                        ,  FunctionTree(*tree.getLeft())); //May be use abstract value instead.
            case FunctionTree::Log:
                break;
            case FunctionTree::Exp:
                return FunctionTree(FunctionTree::Mul
                        , FunctionTree(tree)
                        , derivativeTree(*tree.getLeft()));
            case FunctionTree::Cos: {
                FunctionTree temp(FunctionTree::Mul
                        , FunctionTree(FunctionTree::Sin, FunctionTree(*tree.getLeft()))
                        , derivativeTree(*tree.getLeft()));
                return FunctionTree(FunctionTree::Mul, f.getConstantNode(3), std::move(temp));
            }
            case FunctionTree::Sin:
                return FunctionTree(FunctionTree::Mul
                        , FunctionTree(FunctionTree::Cos, FunctionTree(*tree.getLeft()))
                        , derivativeTree(*tree.getLeft()));
            case FunctionTree::Tan: //TODO
                return FunctionTree(FunctionTree::Mul, derivativeTree(*tree.getLeft()),
                                    FunctionTree(FunctionTree::Square,
                                                 FunctionTree(FunctionTree::Sec, FunctionTree(*tree.getLeft()))));
            case FunctionTree::Sec:
                return FunctionTree(FunctionTree::Mul, derivativeTree(*tree.getLeft()),
                                    FunctionTree(FunctionTree::Mul,
                                                 FunctionTree(FunctionTree::Tan, FunctionTree(*tree.getLeft())),
                                                 FunctionTree(FunctionTree::Sec, FunctionTree(*tree.getLeft()))));
            case FunctionTree::Csc: {
                FunctionTree temp(FunctionTree::Mul,
                             FunctionTree(FunctionTree::Cot, FunctionTree(*tree.getLeft())),
                             FunctionTree(FunctionTree::Csc, FunctionTree(*tree.getLeft())));
                FunctionTree temp1(FunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return FunctionTree(FunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::Cot: {
                FunctionTree temp(FunctionTree::Square
                        , FunctionTree(FunctionTree::Sec, FunctionTree(*tree.getLeft())));
                FunctionTree temp1(FunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return FunctionTree(FunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::ArcCos: {
                FunctionTree square(FunctionTree::Square, FunctionTree(*tree.getLeft()));
                FunctionTree sub(FunctionTree::Sub, f.getConstantNode(1), std::move(square));
                FunctionTree temp(FunctionTree::Reciprocal, FunctionTree(FunctionTree::Sqrt, std::move(sub)));
                FunctionTree temp1(FunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return FunctionTree(FunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::ArcSin: {
                FunctionTree square(FunctionTree::Square, FunctionTree(*tree.getLeft()));
                FunctionTree sub(FunctionTree::Sub, f.getConstantNode(1), std::move(square));
                FunctionTree temp(FunctionTree::Reciprocal, FunctionTree(FunctionTree::Sqrt, std::move(sub)));
                return FunctionTree(FunctionTree::Mul, std::move(temp), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::ArcTan: {
                FunctionTree square(FunctionTree::Square, FunctionTree(*tree.getLeft()));
                FunctionTree add(FunctionTree::Add, f.getConstantNode(1), std::move(square));
                FunctionTree temp(FunctionTree::Reciprocal, std::move(add));
                return FunctionTree(FunctionTree::Mul, std::move(temp), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::ArcSec: {
                FunctionTree square1(FunctionTree::Square, FunctionTree(*tree.getLeft()));
                FunctionTree square2(FunctionTree::Square, FunctionTree(*tree.getLeft()));
                FunctionTree sub(FunctionTree::Sub, std::move(square1), f.getConstantNode(1));
                FunctionTree mul(FunctionTree::Mul, FunctionTree(square2), std::move(sub));
                FunctionTree temp(FunctionTree::Reciprocal, FunctionTree(FunctionTree::Sqrt, std::move(mul)));
                return FunctionTree(FunctionTree::Mul, std::move(temp), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::ArcCsc: {
                FunctionTree square1(FunctionTree::Square, FunctionTree(*tree.getLeft()));
                FunctionTree square2(FunctionTree::Square, FunctionTree(*tree.getLeft()));
                FunctionTree sub(FunctionTree::Sub, std::move(square1), f.getConstantNode(1));
                FunctionTree mul(FunctionTree::Mul, FunctionTree(square2), std::move(sub));
                FunctionTree temp(FunctionTree::Reciprocal, FunctionTree(FunctionTree::Sqrt, std::move(mul)));
                FunctionTree temp1(FunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return FunctionTree(FunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::ArcCot: {
                FunctionTree square(FunctionTree::Square, FunctionTree(*tree.getLeft()));
                FunctionTree add(FunctionTree::Add, f.getConstantNode(1), std::move(square));
                FunctionTree temp(FunctionTree::Reciprocal, std::move(add));
                FunctionTree temp1(FunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return FunctionTree(FunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::Cosh:
                return FunctionTree(FunctionTree::Mul, derivativeTree(*tree.getLeft())
                        , FunctionTree(FunctionTree::Sinh, FunctionTree(*tree.getLeft())));
            case FunctionTree::Sinh:
                return FunctionTree(FunctionTree::Mul, derivativeTree(*tree.getLeft())
                        , FunctionTree(FunctionTree::Cosh, FunctionTree(*tree.getLeft())));
            case FunctionTree::Tanh:
                return FunctionTree(FunctionTree::Mul, derivativeTree(*tree.getLeft())
                        , FunctionTree(FunctionTree::Square
                                , FunctionTree(FunctionTree::Sech, FunctionTree(*tree.getLeft()))));
            case FunctionTree::Sech: {
                FunctionTree temp(FunctionTree::Mul,
                                  FunctionTree(FunctionTree::Sech, FunctionTree(*tree.getLeft())),
                                  FunctionTree(FunctionTree::Tanh, FunctionTree(*tree.getLeft())));
                FunctionTree temp1(FunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return FunctionTree(FunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::Csch: {
                FunctionTree temp(FunctionTree::Mul,
                                  FunctionTree(FunctionTree::Csch, FunctionTree(*tree.getLeft())),
                                  FunctionTree(FunctionTree::Coth, FunctionTree(*tree.getLeft())));
                FunctionTree temp1(FunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return FunctionTree(FunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::Coth: {
                FunctionTree temp(FunctionTree::Square
                        , FunctionTree(FunctionTree::Csch, FunctionTree(*tree.getLeft())));
                FunctionTree temp1(FunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return FunctionTree(FunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case FunctionTree::ArcCosh:
                return FunctionTree(FunctionTree::Div, derivativeTree(*tree.getLeft())
                        , FunctionTree(FunctionTree::Sqrt
                                , FunctionTree(FunctionTree::Sub
                                        , FunctionTree(FunctionTree::Square, FunctionTree(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case FunctionTree::ArcSinh:
                return FunctionTree(FunctionTree::Div, derivativeTree(*tree.getLeft())
                        , FunctionTree(FunctionTree::Sqrt
                                , FunctionTree(FunctionTree::Add
                                        , FunctionTree(FunctionTree::Square, FunctionTree(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case FunctionTree::ArcSech:
                return FunctionTree(FunctionTree::Div
                        , FunctionTree(FunctionTree::Mul, f.getConstantNode(3), derivativeTree(*tree.getLeft()))
                        , FunctionTree(FunctionTree::Mul, FunctionTree(*tree.getLeft())
                                , FunctionTree(FunctionTree::Sub
                                        , FunctionTree(FunctionTree::Square, FunctionTree(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case FunctionTree::ArcCsch:
                return FunctionTree(FunctionTree::Div
                        , FunctionTree(FunctionTree::Mul, f.getConstantNode(3), derivativeTree(*tree.getLeft()))
                        , FunctionTree(FunctionTree::Mul, FunctionTree(*tree.getLeft())
                                , FunctionTree(FunctionTree::Add
                                        , FunctionTree(FunctionTree::Square, FunctionTree(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case FunctionTree::ArcTanh:
            case FunctionTree::ArcCoth:
                return FunctionTree(FunctionTree::Div, derivativeTree(*tree.getLeft())
                        , FunctionTree(FunctionTree::Sub
                                , f.getConstantNode(1)
                                , FunctionTree(FunctionTree::Square, FunctionTree(*tree.getLeft()))));
        }
    }
}