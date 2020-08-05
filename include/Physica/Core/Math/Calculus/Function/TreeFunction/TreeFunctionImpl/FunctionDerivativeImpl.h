/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTIONDERIVATIVEIMPL_H
#define PHYSICA_FUNCTIONDERIVATIVEIMPL_H

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    FunctionDerivative<type, errorTrack>::FunctionDerivative(const TreeFunction<type, errorTrack>& f, size_t index) : f(f), index(index) {}
    /*!
     * Return the partial derivative of the given function.
     */
    template<ScalarType type, bool errorTrack>
    TreeFunction<type, errorTrack> FunctionDerivative<type, errorTrack>::derivative() const {
        const auto variablesLength = f.getVariablesLength();
        const auto constantsLength = f.getConstantsLength() + 4;
        TreeFunction result(derivativeTree(f.getTree()), variablesLength, constantsLength);
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
    template<ScalarType type, bool errorTrack>
    TreeFunctionData<type, errorTrack> FunctionDerivative<type, errorTrack>::derivativeTree(const TreeFunctionData<type, errorTrack>& tree) const {
        const AbstractFunctionTree::FunctionType functionType = tree.getType();
        switch(functionType) {
            case AbstractFunctionTree::Value:
                return index == f.getVariablePos(tree) ? f.getConstantNode(1) : f.getConstantNode(0);
            case AbstractFunctionTree::Add:
            case AbstractFunctionTree::Sub:
                return TreeFunctionData(functionType
                        , derivativeTree(*tree.getLeft()), derivativeTree(*tree.getLeft()));
            case AbstractFunctionTree::Mul:
                return TreeFunctionData(AbstractFunctionTree::Add
                        , TreeFunctionData(functionType, derivativeTree(*tree.getLeft()), TreeFunctionData(*tree.getRight()))
                        , TreeFunctionData(functionType, derivativeTree(*tree.getRight()), TreeFunctionData(*tree.getLeft())));
            case AbstractFunctionTree::Div: {
                TreeFunctionData f1(AbstractFunctionTree::Mul, TreeFunctionData(*tree.getRight()), derivativeTree(*tree.getLeft()));
                TreeFunctionData f2(AbstractFunctionTree::Mul, TreeFunctionData(*tree.getLeft()), derivativeTree(*tree.getRight()));
                TreeFunctionData f1_f2(AbstractFunctionTree::Sub, std::move(f1), std::move(f2));
                TreeFunctionData f2_2(AbstractFunctionTree::Square, TreeFunctionData(*tree.getRight()));
                return TreeFunctionData(AbstractFunctionTree::Div, std::move(f1_f2), std::move(f2_2));
            }
            case AbstractFunctionTree::Square: {
                return TreeFunctionData(AbstractFunctionTree::Mul, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(AbstractFunctionTree::Mul, f.getConstantNode(2), TreeFunctionData(*tree.getLeft())));
            }
            case AbstractFunctionTree::Reciprocal: {
                TreeFunctionData s(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()));
                return TreeFunctionData(AbstractFunctionTree::Mul
                        , f.getConstantNode(3)
                        , TreeFunctionData(AbstractFunctionTree::Div, derivativeTree(*tree.getLeft()), std::move(s)));
            }
            case AbstractFunctionTree::Sqrt: {
                TreeFunctionData s(AbstractFunctionTree::Sqrt, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData r(AbstractFunctionTree::Div, derivativeTree(*tree.getLeft()), std::move(s));
                return TreeFunctionData(AbstractFunctionTree::Div, std::move(r), f.getConstantNode(2));
            }
            case AbstractFunctionTree::Factorial:
                //Error
                return f.getConstantNode(0);
            case AbstractFunctionTree::Ln:
                return TreeFunctionData(AbstractFunctionTree::Div
                        , derivativeTree(*tree.getLeft())
                        , TreeFunctionData(*tree.getLeft())); //May be use abstract value instead.
            case AbstractFunctionTree::Log:
                break;
            case AbstractFunctionTree::Exp:
                return TreeFunctionData(AbstractFunctionTree::Mul
                        , TreeFunctionData(tree)
                        , derivativeTree(*tree.getLeft()));
            case AbstractFunctionTree::Cos: {
                TreeFunctionData temp(AbstractFunctionTree::Mul
                        , TreeFunctionData(AbstractFunctionTree::Sin, TreeFunctionData(*tree.getLeft()))
                        , derivativeTree(*tree.getLeft()));
                return TreeFunctionData(AbstractFunctionTree::Mul, f.getConstantNode(3), std::move(temp));
            }
            case AbstractFunctionTree::Sin:
                return TreeFunctionData(AbstractFunctionTree::Mul
                        , TreeFunctionData(AbstractFunctionTree::Cos, TreeFunctionData(*tree.getLeft()))
                        , derivativeTree(*tree.getLeft()));
            case AbstractFunctionTree::Tan:
                return TreeFunctionData(AbstractFunctionTree::Mul, derivativeTree(*tree.getLeft()),
                                        TreeFunctionData(AbstractFunctionTree::Square,
                                                         TreeFunctionData(AbstractFunctionTree::Sec, TreeFunctionData(*tree.getLeft()))));
            case AbstractFunctionTree::Sec:
                return TreeFunctionData(AbstractFunctionTree::Mul, derivativeTree(*tree.getLeft()),
                                        TreeFunctionData(AbstractFunctionTree::Mul,
                                                         TreeFunctionData(AbstractFunctionTree::Tan, TreeFunctionData(*tree.getLeft())),
                                                         TreeFunctionData(AbstractFunctionTree::Sec, TreeFunctionData(*tree.getLeft()))));
            case AbstractFunctionTree::Csc: {
                TreeFunctionData temp(AbstractFunctionTree::Mul,
                                      TreeFunctionData(AbstractFunctionTree::Cot, TreeFunctionData(*tree.getLeft())),
                                      TreeFunctionData(AbstractFunctionTree::Csc, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(AbstractFunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::Cot: {
                TreeFunctionData temp(AbstractFunctionTree::Square
                        , TreeFunctionData(AbstractFunctionTree::Sec, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(AbstractFunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::ArcCos: {
                TreeFunctionData square(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData sub(AbstractFunctionTree::Sub, f.getConstantNode(1), std::move(square));
                TreeFunctionData temp(AbstractFunctionTree::Reciprocal, TreeFunctionData(AbstractFunctionTree::Sqrt, std::move(sub)));
                TreeFunctionData temp1(AbstractFunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::ArcSin: {
                TreeFunctionData square(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData sub(AbstractFunctionTree::Sub, f.getConstantNode(1), std::move(square));
                TreeFunctionData temp(AbstractFunctionTree::Reciprocal, TreeFunctionData(AbstractFunctionTree::Sqrt, std::move(sub)));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::ArcTan: {
                TreeFunctionData square(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData add(AbstractFunctionTree::Add, f.getConstantNode(1), std::move(square));
                TreeFunctionData temp(AbstractFunctionTree::Reciprocal, std::move(add));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::ArcSec: {
                TreeFunctionData square1(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData square2(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData sub(AbstractFunctionTree::Sub, std::move(square1), f.getConstantNode(1));
                TreeFunctionData mul(AbstractFunctionTree::Mul, TreeFunctionData(square2), std::move(sub));
                TreeFunctionData temp(AbstractFunctionTree::Reciprocal, TreeFunctionData(AbstractFunctionTree::Sqrt, std::move(mul)));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::ArcCsc: {
                TreeFunctionData square1(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData square2(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData sub(AbstractFunctionTree::Sub, std::move(square1), f.getConstantNode(1));
                TreeFunctionData mul(AbstractFunctionTree::Mul, TreeFunctionData(square2), std::move(sub));
                TreeFunctionData temp(AbstractFunctionTree::Reciprocal, TreeFunctionData(AbstractFunctionTree::Sqrt, std::move(mul)));
                TreeFunctionData temp1(AbstractFunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::ArcCot: {
                TreeFunctionData square(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData add(AbstractFunctionTree::Add, f.getConstantNode(1), std::move(square));
                TreeFunctionData temp(AbstractFunctionTree::Reciprocal, std::move(add));
                TreeFunctionData temp1(AbstractFunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::Cosh:
                return TreeFunctionData(AbstractFunctionTree::Mul, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(AbstractFunctionTree::Sinh, TreeFunctionData(*tree.getLeft())));
            case AbstractFunctionTree::Sinh:
                return TreeFunctionData(AbstractFunctionTree::Mul, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(AbstractFunctionTree::Cosh, TreeFunctionData(*tree.getLeft())));
            case AbstractFunctionTree::Tanh:
                return TreeFunctionData(AbstractFunctionTree::Mul, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(AbstractFunctionTree::Square
                                , TreeFunctionData(AbstractFunctionTree::Sech, TreeFunctionData(*tree.getLeft()))));
            case AbstractFunctionTree::Sech: {
                TreeFunctionData temp(AbstractFunctionTree::Mul,
                                      TreeFunctionData(AbstractFunctionTree::Sech, TreeFunctionData(*tree.getLeft())),
                                      TreeFunctionData(AbstractFunctionTree::Tanh, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(AbstractFunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::Csch: {
                TreeFunctionData temp(AbstractFunctionTree::Mul,
                                      TreeFunctionData(AbstractFunctionTree::Csch, TreeFunctionData(*tree.getLeft())),
                                      TreeFunctionData(AbstractFunctionTree::Coth, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(AbstractFunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::Coth: {
                TreeFunctionData temp(AbstractFunctionTree::Square
                        , TreeFunctionData(AbstractFunctionTree::Csch, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(AbstractFunctionTree::Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(AbstractFunctionTree::Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case AbstractFunctionTree::ArcCosh:
                return TreeFunctionData(AbstractFunctionTree::Div, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(AbstractFunctionTree::Sqrt
                                , TreeFunctionData(AbstractFunctionTree::Sub
                                        , TreeFunctionData(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case AbstractFunctionTree::ArcSinh:
                return TreeFunctionData(AbstractFunctionTree::Div, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(AbstractFunctionTree::Sqrt
                                , TreeFunctionData(AbstractFunctionTree::Add
                                        , TreeFunctionData(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case AbstractFunctionTree::ArcSech:
                return TreeFunctionData(AbstractFunctionTree::Div
                        , TreeFunctionData(AbstractFunctionTree::Mul, f.getConstantNode(3), derivativeTree(*tree.getLeft()))
                        , TreeFunctionData(AbstractFunctionTree::Mul, TreeFunctionData(*tree.getLeft())
                                , TreeFunctionData(AbstractFunctionTree::Sub
                                        , TreeFunctionData(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case AbstractFunctionTree::ArcCsch:
                return TreeFunctionData(AbstractFunctionTree::Div
                        , TreeFunctionData(AbstractFunctionTree::Mul, f.getConstantNode(3), derivativeTree(*tree.getLeft()))
                        , TreeFunctionData(AbstractFunctionTree::Mul, TreeFunctionData(*tree.getLeft())
                                , TreeFunctionData(AbstractFunctionTree::Add
                                        , TreeFunctionData(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case AbstractFunctionTree::ArcTanh:
            case AbstractFunctionTree::ArcCoth:
                return TreeFunctionData(AbstractFunctionTree::Div, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(AbstractFunctionTree::Sub
                                , f.getConstantNode(1)
                                , TreeFunctionData(AbstractFunctionTree::Square, TreeFunctionData(*tree.getLeft()))));
        }
    }
}

#endif
