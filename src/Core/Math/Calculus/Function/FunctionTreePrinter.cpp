/*
 * Copyright (c) 2020 NewSigma@163.com. All rights reserved.
 */
#include "Physica/Core/Math/Calculus/Function/FunctionTreePrinter.h"
#include "Physica/Core/Math/Calculus/Function/Function.h"

namespace Physica::Core {
    const char* FunctionTreePrinter::rightBranch = ".---";
    const char* FunctionTreePrinter::leftBranch = "`---";
    const char* FunctionTreePrinter::straightLine = "|";
    const char* FunctionTreePrinter::space = " ";
    const char* FunctionTreePrinter::spaces = "    ";

    FunctionTreePrinter::FunctionTreePrinter(const Function& f, std::ostream& os) : f(f), os(os) {}

    void FunctionTreePrinter::printImpl(const FunctionTree& functionTree, bool isLeft) {
        if(functionTree.type == FunctionTree::Value) {
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
            case FunctionTree::Add:
                os << "+";
                break;
            case FunctionTree::Sub:
                os << "-";
                break;
            case FunctionTree::Mul:
                os << "*";
                break;
            case FunctionTree::Div:
                os << "/";
                break;
            case FunctionTree::Square:
                os << "square";
                break;
            case FunctionTree::Reciprocal:
                os << "reciprocal";
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
                os << "log";
                break;
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

    void FunctionTreePrinter::printList() {
        for(const char* str : list)
            os << str;
    }
}