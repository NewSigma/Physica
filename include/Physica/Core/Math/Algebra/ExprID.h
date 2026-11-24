/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2025 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

namespace Physica {
    /**
     * \class ExprID is used in implementation of expression template and auto differential.
     */
    enum class ExprID : char {
        Minus,
        Add,
        Sub,
        Mul,
        Div,
        MulAdd,
        Sum,
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
        ExpM1,
        Pow,
        Sin,
        Cos,
        Tan,
        Sec,
        ArcCos,
        Cosh,
        Sinh,
        Tanh,
        Sech,
        ArcSinh,
        ArcTanh,
        LnCosh,
        Softmax,
        Sigmoid,
        Softplus,
    };
}
