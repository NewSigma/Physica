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

#include "../Square.h"

namespace Physica {
    template<Vector V>
    void VectorExpr<ExprID::Square, V>::assign_mkl(Vector auto& v) const noexcept {
        using Tm = decltype(std::declval<T>().toMKL());
        v.assert_assign_mkl(Base::getExpr());

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
