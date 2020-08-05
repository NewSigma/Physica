/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_FUNCTIONPRINTERIMPL_H
#define PHYSICA_FUNCTIONPRINTERIMPL_H

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    FunctionPrinter<type, errorTrack>::FunctionPrinter(const TreeFunction<type, errorTrack>& f, std::ostream& os) : f(f), os(os) {}

    template<ScalarType type, bool errorTrack>
    void FunctionPrinter<type, errorTrack>::printImpl(const TreeFunctionData<type, errorTrack>& functionTree) {
        switch(functionTree.getType()) {
            case AbstractFunctionTree::Value: {
                auto pos = f.getVariablePos(functionTree);
                if(pos != 0)
                    os << "x" << pos;
                else
                    os << double(*functionTree.getValue());
                return;
            }
            case AbstractFunctionTree::Add:
                printImpl(*functionTree.getLeft());
                os << " + ";
                printImpl(*functionTree.getRight());
                break;
            case AbstractFunctionTree::Sub:
                printImpl(*functionTree.getLeft());
                os << " - ";
                printImpl(*functionTree.getRight());
                break;
            case AbstractFunctionTree::Mul:
                printImpl(*functionTree.getLeft());
                os << " * ";
                printImpl(*functionTree.getRight());
                break;
            case AbstractFunctionTree::Div:
                printImpl(*functionTree.getLeft());
                os << " / ";
                printImpl(*functionTree.getRight());
                break;
            case AbstractFunctionTree::Square:
                os << '(';
                printImpl(*functionTree.getLeft());
                os << ")^2";
                break;
            case AbstractFunctionTree::Reciprocal:
                os << "1 / ";
                break;
            case AbstractFunctionTree::Sqrt:
                os << "sqrt";
                break;
            case AbstractFunctionTree::Factorial:
                os << "factorial";
                break;
            case AbstractFunctionTree::Ln:
                os << "ln";
                break;
            case AbstractFunctionTree::Log:
                os << "log_(";
                printImpl(*functionTree.getLeft());
                os << ")^(";
                printImpl(*functionTree.getRight());
                os << ')';
                return;
            case AbstractFunctionTree::Exp:
                os << "exp";
                break;
            case AbstractFunctionTree::Cos:
                os << "cos";
                break;
            case AbstractFunctionTree::Sin:
                os << "sin";
                break;
            case AbstractFunctionTree::Tan:
                os << "tan";
                break;
            case AbstractFunctionTree::Sec:
                os << "sec";
                break;
            case AbstractFunctionTree::Csc:
                os << "csc";
                break;
            case AbstractFunctionTree::Cot:
                os << "cot";
                break;
            case AbstractFunctionTree::ArcCos:
                os << "arccos";
                break;
            case AbstractFunctionTree::ArcSin:
                os << "arcsin";
                break;
            case AbstractFunctionTree::ArcTan:
                os << "arctan";
                break;
            case AbstractFunctionTree::ArcSec:
                os << "arcsec";
                break;
            case AbstractFunctionTree::ArcCsc:
                os << "arccsc";
                break;
            case AbstractFunctionTree::ArcCot:
                os << "arccot";
                break;
            case AbstractFunctionTree::Cosh:
                os << "cosh";
                break;
            case AbstractFunctionTree::Sinh:
                os << "sinh";
                break;
            case AbstractFunctionTree::Tanh:
                os << "tanh";
                break;
            case AbstractFunctionTree::Sech:
                os << "sech";
                break;
            case AbstractFunctionTree::Csch:
                os << "csch";
                break;
            case AbstractFunctionTree::Coth:
                os << "coth";
                break;
            case AbstractFunctionTree::ArcCosh:
                os << "arccosh";
                break;
            case AbstractFunctionTree::ArcSinh:
                os << "arcsinh";
                break;
            case AbstractFunctionTree::ArcTanh:
                os << "arctanh";
                break;
            case AbstractFunctionTree::ArcSech:
                os << "arcsech";
                break;
            case AbstractFunctionTree::ArcCsch:
                os << "arccsch";
                break;
            case AbstractFunctionTree::ArcCoth:
                os << "arccoth";
                break;
        }
        os << '(';
        printImpl(*functionTree.getLeft());
        os << ')';
    }
}

#endif
