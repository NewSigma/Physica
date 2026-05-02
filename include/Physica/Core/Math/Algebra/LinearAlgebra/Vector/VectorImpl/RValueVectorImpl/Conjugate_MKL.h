/*
 * Copyright 2025-2026 Weibo He.
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
