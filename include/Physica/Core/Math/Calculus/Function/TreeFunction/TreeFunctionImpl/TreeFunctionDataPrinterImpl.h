/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#ifndef PHYSICA_TREEFUNCTIONDATAPRINTERIMPL_H
#define PHYSICA_TREEFUNCTIONDATAPRINTERIMPL_H

namespace Physica::Core {
    template<ScalarType type, bool errorTrack>
    TreeFunctionPrinter<type, errorTrack>::TreeFunctionPrinter(const TreeFunction<type, errorTrack>& f, std::ostream& os) : f(f), os(os) {}

    template<ScalarType type, bool errorTrack>
    void TreeFunctionPrinter<type, errorTrack>::printImpl(const TreeFunctionData<type, errorTrack>& functionTree, bool isLeft) {
        if(functionTree.getType() == Value) {
            auto pos = f.getVariablePos(functionTree);
            if(pos != 0) {
                printList();
                os << (isLeft ? leftBranch : rightBranch) << "x" << pos << '\n';
            }
            else
                os << double(*functionTree.value);
            return;
        }
        const bool isBranch = &functionTree != &f.getTree();
        const bool rightExist = functionTree.right != nullptr;
        //Right
        if(isBranch)
            list.push_back(isLeft ? straightLine : space);
        list.push_back(spaces);

        if(rightExist)
            printImpl(*functionTree.right, false);

        if(isBranch)
            list.pop_back();
        list.pop_back();
        //Center
        list.push_back(space);
        list.push_back(isBranch
                       ? (isLeft ? leftBranch : rightBranch)
                       : spaces);
        printList();
        list.pop_back();
        list.pop_back();
        switch(functionTree.type) {
            case Add:
                os << "+";
                break;
            case Sub:
                os << "-";
                break;
            case Mul:
                os << "*";
                break;
            case Div:
                os << "/";
                break;
            case Square:
                os << "square";
                break;
            case Reciprocal:
                os << "reciprocal";
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
                os << "log";
                break;
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
            default:;
        }
        os << '\n';
        //Left
        if(isBranch)
            list.push_back(isLeft ? space : straightLine);
        list.push_back(spaces);
        printImpl(*functionTree.left, true);
        if(isBranch)
            list.pop_back();
        list.pop_back();
    }

    template<ScalarType type, bool errorTrack>
    void TreeFunctionPrinter<type, errorTrack>::printList() {
        for(const char* str : list)
            os << str;
    }
}

#endif
