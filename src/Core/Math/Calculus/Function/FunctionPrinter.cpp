/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/Math/Calculus/Function/FunctionPrinter.h"

namespace Physica::Core {
    FunctionPrinter::FunctionPrinter(const Function& f, std::ostream& os) : f(f), os(os) {}

    void FunctionPrinter::printImpl(const FunctionTree& functionTree) {
        switch(functionTree.getType()) {
            case FunctionTree::Value: {
                auto pos = f.getVariablePos(functionTree);
                if(pos != 0)
                    os << "x" << pos;
                else
                    os << double(*functionTree.getValue());
                return;
            }
            case FunctionTree::Add:
                printImpl(*functionTree.getLeft());
                os << " + ";
                printImpl(*functionTree.getRight());
                break;
            case FunctionTree::Sub:
                printImpl(*functionTree.getLeft());
                os << " - ";
                printImpl(*functionTree.getRight());
                break;
            case FunctionTree::Mul:
                printImpl(*functionTree.getLeft());
                os << " * ";
                printImpl(*functionTree.getRight());
                break;
            case FunctionTree::Div:
                printImpl(*functionTree.getLeft());
                os << " / ";
                printImpl(*functionTree.getRight());
                break;
            case FunctionTree::Square:
                os << '(';
                printImpl(*functionTree.getLeft());
                os << ")^2";
                break;
            case FunctionTree::Reciprocal:
                os << "1 / ";
                break;
            case FunctionTree::Sqrt:
                os << "sqrt";
                break;
            case FunctionTree::Factorial:
                os << "factorial";
                break;
            case FunctionTree::Ln:
                os << "ln";
                break;
            case FunctionTree::Log:
                os << "log_(";
                printImpl(*functionTree.getLeft());
                os << ")^(";
                printImpl(*functionTree.getRight());
                os << ')';
                return;
            case FunctionTree::Exp:
                os << "exp";
                break;
            case FunctionTree::Cos:
                os << "cos";
                break;
            case FunctionTree::Sin:
                os << "sin";
                break;
            case FunctionTree::Tan:
                os << "tan";
                break;
            case FunctionTree::Sec:
                os << "sec";
                break;
            case FunctionTree::Csc:
                os << "csc";
                break;
            case FunctionTree::Cot:
                os << "cot";
                break;
            case FunctionTree::ArcCos:
                os << "arccos";
                break;
            case FunctionTree::ArcSin:
                os << "arcsin";
                break;
            case FunctionTree::ArcTan:
                os << "arctan";
                break;
            case FunctionTree::ArcSec:
                os << "arcsec";
                break;
            case FunctionTree::ArcCsc:
                os << "arccsc";
                break;
            case FunctionTree::ArcCot:
                os << "arccot";
                break;
            case FunctionTree::Cosh:
                os << "cosh";
                break;
            case FunctionTree::Sinh:
                os << "sinh";
                break;
            case FunctionTree::Tanh:
                os << "tanh";
                break;
            case FunctionTree::Sech:
                os << "sech";
                break;
            case FunctionTree::Csch:
                os << "csch";
                break;
            case FunctionTree::Coth:
                os << "coth";
                break;
            case FunctionTree::ArcCosh:
                os << "arccosh";
                break;
            case FunctionTree::ArcSinh:
                os << "arcsinh";
                break;
            case FunctionTree::ArcTanh:
                os << "arctanh";
                break;
            case FunctionTree::ArcSech:
                os << "arcsech";
                break;
            case FunctionTree::ArcCsch:
                os << "arccsch";
                break;
            case FunctionTree::ArcCoth:
                os << "arccoth";
                break;
        }
        os << '(';
        printImpl(*functionTree.getLeft());
        os << ')';
    }
}