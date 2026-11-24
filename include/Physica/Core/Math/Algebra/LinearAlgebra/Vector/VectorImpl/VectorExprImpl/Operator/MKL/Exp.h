/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include "../Exp.h"

namespace Physica {
    template<Vector V>
    void VectorExpr<ExprID::Exp, V>::assign_mkl(Vector auto& v) const noexcept {
        using Tm = decltype(std::declval<T>().toMKL());
        v.assert_assign_lapack(Base::getExpr());

        size_t n = Base::getLength();
        const auto* a = reinterpret_cast<const Tm*>(Base::getExpr().data_handle());
        auto* y = reinterpret_cast<Tm*>(v.data_handle());
        if constexpr (isComplex()) {
            if constexpr (T::Prec == Float32)
                vcExp_64(n, a, y);
            else
                vzExp_64(n, a, y);
        }
        else {
            if constexpr (T::Prec == Float16)
                vhExp_64(n, a, y);
            else if constexpr (T::Prec == Float32)
                vsExp_64(n, a, y);
            else
                vdExp_64(n, a, y);
        }
    }
}
