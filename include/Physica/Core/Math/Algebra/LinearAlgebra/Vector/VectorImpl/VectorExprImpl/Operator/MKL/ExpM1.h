/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../ExpM1.h"

namespace Physica {
    template<Vector V>
    void VectorExpr<ExprID::ExpM1, V>::assign_mkl(Vector auto& v) const noexcept {
        static_assert(!isComplex(), "[Error]: MKL has no complex expm1");
        using Tm = decltype(std::declval<T>().toMKL());
        v.assert_assign_lapack(Base::getExpr());

        size_t n = Base::getLength();
        const auto* a = reinterpret_cast<const Tm*>(Base::getExpr().data());
        auto* y = reinterpret_cast<Tm*>(v.data());
        if constexpr (T::Prec == Float32)
            vsExpm1_64(n, a, y);
        else
            vdExpm1_64(n, a, y);
    }
}
