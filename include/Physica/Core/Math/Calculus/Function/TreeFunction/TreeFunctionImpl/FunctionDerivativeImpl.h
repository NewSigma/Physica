/*
 * Copyright 2020 WeiBo He.
 *
 * This file is part of Physica.

 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef PHYSICA_FUNCTIONDERIVATIVEIMPL_H
#define PHYSICA_FUNCTIONDERIVATIVEIMPL_H

namespace Physica::Core {
    template<ScalarOption option, bool errorTrack>
    FunctionDerivative<option, errorTrack>::FunctionDerivative(const TreeFunctionData<option, errorTrack>& f, size_t index) : f(f), index(index) {}
    /*!
     * Return the partial derivative of the given function.
     */
    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack> FunctionDerivative<option, errorTrack>::derivative() const {
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
    template<ScalarOption option, bool errorTrack>
    TreeFunctionData<option, errorTrack> FunctionDerivative<option, errorTrack>::derivativeTree(const TreeFunctionData<option, errorTrack>& tree) const {
        const FunctionType functionType = tree.getType();
        switch(functionType) {
            case Value:
                return index == f.getVariablePos(tree) ? f.getConstantNode(1) : f.getConstantNode(0);
            case Add:
            case Sub:
                return TreeFunctionData(functionType
                        , derivativeTree(*tree.getLeft()), derivativeTree(*tree.getLeft()));
            case Mul:
                return TreeFunctionData(Add
                        , TreeFunctionData(functionType, derivativeTree(*tree.getLeft()), TreeFunctionData(*tree.getRight()))
                        , TreeFunctionData(functionType, derivativeTree(*tree.getRight()), TreeFunctionData(*tree.getLeft())));
            case Div: {
                TreeFunctionData f1(Mul, TreeFunctionData(*tree.getRight()), derivativeTree(*tree.getLeft()));
                TreeFunctionData f2(Mul, TreeFunctionData(*tree.getLeft()), derivativeTree(*tree.getRight()));
                TreeFunctionData f1_f2(Sub, std::move(f1), std::move(f2));
                TreeFunctionData f2_2(Square, TreeFunctionData(*tree.getRight()));
                return TreeFunctionData(Div, std::move(f1_f2), std::move(f2_2));
            }
            case Square: {
                return TreeFunctionData(Mul, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(Mul, f.getConstantNode(2), TreeFunctionData(*tree.getLeft())));
            }
            case Reciprocal: {
                TreeFunctionData s(Square, TreeFunctionData(*tree.getLeft()));
                return TreeFunctionData(Mul
                        , f.getConstantNode(3)
                        , TreeFunctionData(Div, derivativeTree(*tree.getLeft()), std::move(s)));
            }
            case Sqrt: {
                TreeFunctionData s(Sqrt, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData r(Div, derivativeTree(*tree.getLeft()), std::move(s));
                return TreeFunctionData(Div, std::move(r), f.getConstantNode(2));
            }
            case Factorial:
                //Error
                return f.getConstantNode(0);
            case Ln:
                return TreeFunctionData(Div
                        , derivativeTree(*tree.getLeft())
                        , TreeFunctionData(*tree.getLeft())); //May be use abstract value instead.
            case Log:
                break;
            case Exp:
                return TreeFunctionData(Mul
                        , TreeFunctionData(tree)
                        , derivativeTree(*tree.getLeft()));
            case Cos: {
                TreeFunctionData temp(Mul
                        , TreeFunctionData(Sin, TreeFunctionData(*tree.getLeft()))
                        , derivativeTree(*tree.getLeft()));
                return TreeFunctionData(Mul, f.getConstantNode(3), std::move(temp));
            }
            case Sin:
                return TreeFunctionData(Mul
                        , TreeFunctionData(Cos, TreeFunctionData(*tree.getLeft()))
                        , derivativeTree(*tree.getLeft()));
            case Tan:
                return TreeFunctionData(Mul, derivativeTree(*tree.getLeft()),
                                        TreeFunctionData(Square,
                                                         TreeFunctionData(Sec, TreeFunctionData(*tree.getLeft()))));
            case Sec:
                return TreeFunctionData(Mul, derivativeTree(*tree.getLeft()),
                                        TreeFunctionData(Mul,
                                                         TreeFunctionData(Tan, TreeFunctionData(*tree.getLeft())),
                                                         TreeFunctionData(Sec, TreeFunctionData(*tree.getLeft()))));
            case Csc: {
                TreeFunctionData temp(Mul,
                                      TreeFunctionData(Cot, TreeFunctionData(*tree.getLeft())),
                                      TreeFunctionData(Csc, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case Cot: {
                TreeFunctionData temp(Square
                        , TreeFunctionData(Sec, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case ArcCos: {
                TreeFunctionData square(Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData sub(Sub, f.getConstantNode(1), std::move(square));
                TreeFunctionData temp(Reciprocal, TreeFunctionData(Sqrt, std::move(sub)));
                TreeFunctionData temp1(Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case ArcSin: {
                TreeFunctionData square(Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData sub(Sub, f.getConstantNode(1), std::move(square));
                TreeFunctionData temp(Reciprocal, TreeFunctionData(Sqrt, std::move(sub)));
                return TreeFunctionData(Mul, std::move(temp), derivativeTree(*tree.getLeft()));
            }
            case ArcTan: {
                TreeFunctionData square(Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData add(Add, f.getConstantNode(1), std::move(square));
                TreeFunctionData temp(Reciprocal, std::move(add));
                return TreeFunctionData(Mul, std::move(temp), derivativeTree(*tree.getLeft()));
            }
            case ArcSec: {
                TreeFunctionData square1(Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData square2(Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData sub(Sub, std::move(square1), f.getConstantNode(1));
                TreeFunctionData mul(Mul, TreeFunctionData(square2), std::move(sub));
                TreeFunctionData temp(Reciprocal, TreeFunctionData(Sqrt, std::move(mul)));
                return TreeFunctionData(Mul, std::move(temp), derivativeTree(*tree.getLeft()));
            }
            case ArcCsc: {
                TreeFunctionData square1(Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData square2(Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData sub(Sub, std::move(square1), f.getConstantNode(1));
                TreeFunctionData mul(Mul, TreeFunctionData(square2), std::move(sub));
                TreeFunctionData temp(Reciprocal, TreeFunctionData(Sqrt, std::move(mul)));
                TreeFunctionData temp1(Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case ArcCot: {
                TreeFunctionData square(Square, TreeFunctionData(*tree.getLeft()));
                TreeFunctionData add(Add, f.getConstantNode(1), std::move(square));
                TreeFunctionData temp(Reciprocal, std::move(add));
                TreeFunctionData temp1(Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case Cosh:
                return TreeFunctionData(Mul, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(Sinh, TreeFunctionData(*tree.getLeft())));
            case Sinh:
                return TreeFunctionData(Mul, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(Cosh, TreeFunctionData(*tree.getLeft())));
            case Tanh:
                return TreeFunctionData(Mul, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(Square
                                , TreeFunctionData(Sech, TreeFunctionData(*tree.getLeft()))));
            case Sech: {
                TreeFunctionData temp(Mul,
                                      TreeFunctionData(Sech, TreeFunctionData(*tree.getLeft())),
                                      TreeFunctionData(Tanh, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case Csch: {
                TreeFunctionData temp(Mul,
                                      TreeFunctionData(Csch, TreeFunctionData(*tree.getLeft())),
                                      TreeFunctionData(Coth, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case Coth: {
                TreeFunctionData temp(Square
                        , TreeFunctionData(Csch, TreeFunctionData(*tree.getLeft())));
                TreeFunctionData temp1(Mul, f.getConstantNode(3), std::move(temp));
                return TreeFunctionData(Mul, std::move(temp1), derivativeTree(*tree.getLeft()));
            }
            case ArcCosh:
                return TreeFunctionData(Div, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(Sqrt
                                , TreeFunctionData(Sub
                                        , TreeFunctionData(Square, TreeFunctionData(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case ArcSinh:
                return TreeFunctionData(Div, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(Sqrt
                                , TreeFunctionData(Add
                                        , TreeFunctionData(Square, TreeFunctionData(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case ArcSech:
                return TreeFunctionData(Div
                        , TreeFunctionData(Mul, f.getConstantNode(3), derivativeTree(*tree.getLeft()))
                        , TreeFunctionData(Mul, TreeFunctionData(*tree.getLeft())
                                , TreeFunctionData(Sub
                                        , TreeFunctionData(Square, TreeFunctionData(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case ArcCsch:
                return TreeFunctionData(Div
                        , TreeFunctionData(Mul, f.getConstantNode(3), derivativeTree(*tree.getLeft()))
                        , TreeFunctionData(Mul, TreeFunctionData(*tree.getLeft())
                                , TreeFunctionData(Add
                                        , TreeFunctionData(Square, TreeFunctionData(*tree.getLeft()))
                                        , f.getConstantNode(1))));
            case ArcTanh:
            case ArcCoth:
                return TreeFunctionData(Div, derivativeTree(*tree.getLeft())
                        , TreeFunctionData(Sub
                                , f.getConstantNode(1)
                                , TreeFunctionData(Square, TreeFunctionData(*tree.getLeft()))));
        }
    }
}

#endif
