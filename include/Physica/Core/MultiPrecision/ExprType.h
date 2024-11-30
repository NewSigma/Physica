/*
 * Copyright 2021-2024 Weibo He.
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
#pragma once

namespace Physica::Core {
    /**
     * \class ExprType is used in implementation of expression template and auto differential.
     */
    enum class ExprType : char {
        Set,
        Assign, //Cherry pick a node at forward pass, collect grads at backward pass
        Diff, //Declare a differential at forward pass, jump to the grad node at backward pass
        Minus,
        Add,
        Sub,
        Mul,
        Div,
        Sum,
        MulAdd2,
        MulAdd4,
        MulAdd8,
        More,
        MoreEq,
        Reciprocal,
        Sqrt,
        Cbrt,
        Abs,
        Relu,
        Unit,
        Square,
        Ln,
        Ln1p,
        Exp,
        Pow,
        Sin,
        Cos,
        Tan,
        Sec,
        ArcCos,
        Cosh,
        Tanh,
        Sech,
        LnCosh
    };

    constexpr const char* exprTypeToStr(ExprType type) {
        using enum ExprType;
        switch (type) {
            case Set: return "Set";
            case Assign: return "Assign";
            case Diff: return "Diff";
            case Minus: return "Minus";
            case Add: return "Add";
            case Sub: return "Sub";
            case Mul: return "Mul";
            case Div: return "Div";
            case Sum: return "Sum";
            case MulAdd2: return "MulAdd2";
            case MulAdd4: return "MulAdd4";
            case MulAdd8: return "MulAdd8";
            case More: return "More";
            case MoreEq: return "MoreEq";
            case Reciprocal: return "Reciprocal";
            case Sqrt: return "Sqrt";
            case Cbrt: return "Cbrt";
            case Abs: return "Abs";
            case Relu: return "Relu";
            case Unit: return "Unit";
            case Square: return "Square";
            case Ln: return "Ln";
            case Ln1p: return "Ln1p";
            case Exp: return "Exp";
            case Pow: return "Pow";
            case Sin: return "Sin";
            case Cos: return "Cos";
            case Tan: return "Tan";
            case Sec: return "Sec";
            case ArcCos: return "ArcCos";
            case Cosh: return "Cosh";
            case Tanh: return "Tanh";
            case Sech: return "Sech";
            case LnCosh: return "LnCosh";
            default: [[unlikely]]
                return "Unknown";
        }
    }
}
