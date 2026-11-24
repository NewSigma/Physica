/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Physica/Core/Math/Algebra/ExprID.h"

namespace Physica::Internal {
    /**
     * \class UnaryScalarOpRtnTy returns the scalar type of a unary expression
     * identified by \tparam ID applied to scalar type \tparam T.
     */
    template<ExprID ID, class T>
    struct UnaryScalarOpRtnTy {
        using Type = T;
    };

    template<class T>
    struct UnaryScalarOpRtnTy<ExprID::Abs, T> {
        using Type = T::RealType;
    };

    template<class T>
    struct UnaryScalarOpRtnTy<ExprID::Unit, T> {
        using Type = T::ValueType;
    };
}
