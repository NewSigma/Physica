/*
 * Copyright 2021-2024 WeiBo He.
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
     * \class ExpressionType is used in implementation of expression template and auto differential.
     */
    enum class ExpressionType : char {
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
        Exp,
        Pow,
        Sin,
        Cos,
        ArcCos,
        Cosh
    };

    constexpr const char* expressionTypeToStr(ExpressionType type) {
        switch (type) {
            case ExpressionType::Set: return "Set";
            case ExpressionType::Assign: return "Assign";
            case ExpressionType::Diff: return "Diff";
            case ExpressionType::Minus: return "Minus";
            case ExpressionType::Add: return "Add";
            case ExpressionType::Sub: return "Sub";
            case ExpressionType::Mul: return "Mul";
            case ExpressionType::Div: return "Div";
            case ExpressionType::Sum: return "Sum";
            case ExpressionType::MulAdd2: return "MulAdd2";
            case ExpressionType::MulAdd4: return "MulAdd4";
            case ExpressionType::MulAdd8: return "MulAdd8";
            case ExpressionType::More: return "More";
            case ExpressionType::MoreEq: return "MoreEq";
            case ExpressionType::Reciprocal: return "Reciprocal";
            case ExpressionType::Sqrt: return "Sqrt";
            case ExpressionType::Cbrt: return "Cbrt";
            case ExpressionType::Abs: return "Abs";
            case ExpressionType::Relu: return "Relu";
            case ExpressionType::Square: return "Square";
            case ExpressionType::Ln: return "Ln";
            case ExpressionType::Exp: return "Exp";
            case ExpressionType::Pow: return "Pow";
            case ExpressionType::Sin: return "Sin";
            case ExpressionType::Cos: return "Cos";
            case ExpressionType::ArcCos: return "ArcCos";
            default: [[unlikely]]
                return "Unknown";
        }
    }
}
