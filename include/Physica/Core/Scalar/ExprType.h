/*
 * Copyright 2021-2025 Weibo He.
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

namespace Physica {
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
        MulAdd,
        Sum,
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
        Ln1pExp,
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
        ArcSinh,
        ArcTanh,
        LnCosh,
        Softmax
    };
}
