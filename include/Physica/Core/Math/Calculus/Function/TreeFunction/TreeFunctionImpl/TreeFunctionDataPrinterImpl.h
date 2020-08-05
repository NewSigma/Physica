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
        if(functionTree.getType() == AbstractFunctionTree::Value) {
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
            case AbstractFunctionTree::Add:
                os << "+";
                break;
            case AbstractFunctionTree::Sub:
                os << "-";
                break;
            case AbstractFunctionTree::Mul:
                os << "*";
                break;
            case AbstractFunctionTree::Div:
                os << "/";
                break;
            case AbstractFunctionTree::Square:
                os << "square";
                break;
            case AbstractFunctionTree::Reciprocal:
                os << "reciprocal";
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
                os << "log";
                break;
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
