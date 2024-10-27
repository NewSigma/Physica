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
        ArcCos,
        Cosh,
        Tanh,
        Sech,
        LnCosh
    };

    constexpr const char* exprTypeToStr(ExprType type) {
        switch (type) {
            case ExprType::Set: return "Set";
            case ExprType::Assign: return "Assign";
            case ExprType::Diff: return "Diff";
            case ExprType::Minus: return "Minus";
            case ExprType::Add: return "Add";
            case ExprType::Sub: return "Sub";
            case ExprType::Mul: return "Mul";
            case ExprType::Div: return "Div";
            case ExprType::Sum: return "Sum";
            case ExprType::MulAdd2: return "MulAdd2";
            case ExprType::MulAdd4: return "MulAdd4";
            case ExprType::MulAdd8: return "MulAdd8";
            case ExprType::More: return "More";
            case ExprType::MoreEq: return "MoreEq";
            case ExprType::Reciprocal: return "Reciprocal";
            case ExprType::Sqrt: return "Sqrt";
            case ExprType::Cbrt: return "Cbrt";
            case ExprType::Abs: return "Abs";
            case ExprType::Relu: return "Relu";
            case ExprType::Unit: return "Unit";
            case ExprType::Square: return "Square";
            case ExprType::Ln: return "Ln";
            case ExprType::Ln1p: return "Ln1p";
            case ExprType::Exp: return "Exp";
            case ExprType::Pow: return "Pow";
            case ExprType::Sin: return "Sin";
            case ExprType::Cos: return "Cos";
            case ExprType::ArcCos: return "ArcCos";
            case ExprType::Cosh: return "Cosh";
            case ExprType::Tanh: return "Tanh";
            case ExprType::Sech: return "Sech";
            case ExprType::LnCosh: return "LnCosh";
            default: [[unlikely]]
                return "Unknown";
        }
    }
}
