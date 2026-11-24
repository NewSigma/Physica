/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../Square.h"

namespace Physica {
    template<Vector V>
    void VectorExpr<ExprID::Square, V>::assign_mkl(Vector auto& v) const noexcept {
        using Tm = decltype(std::declval<T>().toMKL());
        v.assert_assign_lapack(Base::getExpr());

        size_t n = Base::getLength();
        const auto* a = reinterpret_cast<const Tm*>(Base::getExpr().data());
        auto* y = reinterpret_cast<Tm*>(v.data());

        if constexpr (T::Prec == Float16)
            vhSqr_64(n, a, y);
        else if constexpr (T::Prec == Float32)
            vsSqr_64(n, a, y);
        else
            vdSqr_64(n, a, y);
    }
}
