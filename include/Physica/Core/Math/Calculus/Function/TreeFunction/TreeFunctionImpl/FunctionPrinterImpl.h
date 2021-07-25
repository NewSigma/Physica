/*
 * Copyright 2020-2021 WeiBo He.
 *
 * This file is part of Physica.
 *
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
#ifndef PHYSICA_FUNCTIONPRINTERIMPL_H
#define PHYSICA_FUNCTIONPRINTERIMPL_H

namespace Physica::Core {
    template<ScalarOption option, bool errorTrack>
    FunctionPrinter<option, errorTrack>::FunctionPrinter(const TreeFunctionData<option, errorTrack>& f_, std::ostream& os) : f(f_), os(os) {}

    template<ScalarOption option, bool errorTrack>
    void FunctionPrinter<option, errorTrack>::printImpl(const TreeFunctionData<option, errorTrack>& functionTree) {
        switch(functionTree.getType()) {
            case Value: {
                auto pos = f.getVariablePos(functionTree);
                if(pos != 0)
                    os << "x" << pos;
                else
                    os << double(*functionTree.getValue(f));
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
