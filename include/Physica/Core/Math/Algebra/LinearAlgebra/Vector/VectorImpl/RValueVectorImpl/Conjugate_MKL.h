/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2025-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

#include "Conjugate.h"

namespace Physica {
    template<Vector V>
    void Conjugate<V>::assign_mkl(Vector auto& v) const noexcept {
        using Tm = decltype(std::declval<T>().toMKL());
        v.assert_assign_lapack(vec);

        size_t n = Base::getLength();
        const auto* a = reinterpret_cast<const Tm*>(vec.data());
        auto* y = reinterpret_cast<Tm*>(v.data());
        if constexpr (T::Prec == Float32)
            vcConj_64(n, a, y);
        else
            vzConj_64(n, a, y);
    }
}
