/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "../Add.h"

namespace Physica {
    template<Vector V1, Vector V2>
    void VectorExpr<ExprID::Add, V1, V2>::assign_mkl(Vector auto& v) const noexcept {
        using Tm = decltype(std::declval<T>().toMKL());
        v.assert_assign_lapack(Base::getLHS());
        v.assert_assign_lapack(Base::getRHS());

        size_t n = Base::getLength();
        const auto* a = reinterpret_cast<const Tm*>(Base::getLHS().data_handle());
        const auto* b = reinterpret_cast<const Tm*>(Base::getRHS().data_handle());
        auto* r = reinterpret_cast<Tm*>(v.data_handle());
        if constexpr (isComplex()) {
            if constexpr (T::Prec == Float32)
                vcAdd_64(n, a, b, r);
            else
                vzAdd_64(n, a, b, r);
        }
        else {
            if constexpr (T::Prec == Float16)
                vhAdd_64(n, a, b, r);
            else if constexpr (T::Prec == Float32)
                vsAdd_64(n, a, b, r);
            else
                vdAdd_64(n, a, b, r);
        }
    }
}
