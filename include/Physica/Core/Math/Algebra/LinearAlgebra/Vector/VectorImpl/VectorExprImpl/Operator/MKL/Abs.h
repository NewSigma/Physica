/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../Abs.h"

namespace Physica {
    template<Vector V>
    void VectorExpr<ExprID::Abs, V>::assign_mkl(Vector auto& v) const noexcept {
        using Tm = decltype(std::declval<T>().toMKL());
        v.assert_assign_lapack(Base::getExpr());

        size_t n = Base::getLength();
        const auto* a = reinterpret_cast<const Tm*>(Base::getExpr().data());
        auto* y = reinterpret_cast<Tm*>(v.data());
        if constexpr (isComplexV) {
            if constexpr (T::Prec == Float32)
                vcAbs_64(n, a, y);
            else
                vzAbs_64(n, a, y);
        }
        else {
            if constexpr (T::Prec == Float16)
                vhAbs_64(n, a, y);
            else if constexpr (T::Prec == Float32)
                vsAbs_64(n, a, y);
            else
                vdAbs_64(n, a, y);
        }
    }
}
