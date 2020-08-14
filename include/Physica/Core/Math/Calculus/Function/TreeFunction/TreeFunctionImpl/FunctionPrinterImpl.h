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
            case Value: {
                auto pos = f.getVariablePos(functionTree);
                if(pos != 0)
                    os << "x" << pos;
                else
                    os << double(*functionTree.getValue());
                return;
            }
            case Add:
                printImpl(*functionTree.getLeft());
                os << " + ";
                printImpl(*functionTree.getRight());
                break;
            case Sub:
                printImpl(*functionTree.getLeft());
                os << " - ";
                printImpl(*functionTree.getRight());
                break;
            case Mul:
                printImpl(*functionTree.getLeft());
                os << " * ";
                printImpl(*functionTree.getRight());
                break;
            case Div:
                printImpl(*functionTree.getLeft());
                os << " / ";
                printImpl(*functionTree.getRight());
                break;
            case Square:
                os << '(';
                printImpl(*functionTree.getLeft());
                os << ")^2";
                break;
            case Reciprocal:
                os << "1 / ";
                break;
            case Sqrt:
                os << "sqrt";
                break;
            case Factorial:
                os << "factorial";
                break;
            case Ln:
                os << "ln";
                break;
            case Log:
                os << "log_(";
                printImpl(*functionTree.getLeft());
                os << ")^(";
                printImpl(*functionTree.getRight());
                os << ')';
                return;
            case Exp:
                os << "exp";
                break;
            case Cos:
                os << "cos";
                break;
            case Sin:
                os << "sin";
                break;
            case Tan:
                os << "tan";
                break;
            case Sec:
                os << "sec";
                break;
            case Csc:
                os << "csc";
                break;
            case Cot:
                os << "cot";
                break;
            case ArcCos:
                os << "arccos";
                break;
            case ArcSin:
                os << "arcsin";
                break;
            case ArcTan:
                os << "arctan";
                break;
            case ArcSec:
                os << "arcsec";
                break;
            case ArcCsc:
                os << "arccsc";
                break;
            case ArcCot:
                os << "arccot";
                break;
            case Cosh:
                os << "cosh";
                break;
            case Sinh:
                os << "sinh";
                break;
            case Tanh:
                os << "tanh";
                break;
            case Sech:
                os << "sech";
                break;
            case Csch:
                os << "csch";
                break;
            case Coth:
                os << "coth";
                break;
            case ArcCosh:
                os << "arccosh";
                break;
            case ArcSinh:
                os << "arcsinh";
                break;
            case ArcTanh:
                os << "arctanh";
                break;
            case ArcSech:
                os << "arcsech";
                break;
            case ArcCsch:
                os << "arccsch";
                break;
            case ArcCoth:
                os << "arccoth";
                break;
        }
        os << '(';
        printImpl(*functionTree.getLeft());
        os << ')';
    }
}

#endif
