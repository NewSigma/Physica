/*
 * Copyright 2026 Weibo He.
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
