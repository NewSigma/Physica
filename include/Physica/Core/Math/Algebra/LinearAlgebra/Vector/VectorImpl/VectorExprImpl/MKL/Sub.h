/*
 * Copyright 2025 Weibo He.
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

#include "../Sub.h"

namespace Physica {
    template<Vector V1, Vector V2>
    void VectorExpr<ExprID::Sub, V1, V2>::assign_mkl(Vector auto& v) const noexcept {
        using Tm = std::conditional<isComplex, typename Tc::MKL_Complex, typename T::MachineType>::type;
        v.assert_assign_mkl(Base::getLHS());
        v.assert_assign_mkl(Base::getRHS());

        size_t n = Base::getLength();
        const auto* a = reinterpret_cast<const Tm*>(Base::getLHS().data());
        const auto* b = reinterpret_cast<const Tm*>(Base::getRHS().data());
        auto* r = reinterpret_cast<Tm*>(v.data());
        if constexpr (isComplex) {
            if constexpr (T::Prec == Float32)
                vcSub_64(n, a, b, r);
            else
                vzSub_64(n, a, b, r);
        }
        else {
            if constexpr (T::Prec == Float16)
                vhSub_64(n, a, b, r);
            else if constexpr (T::Prec == Float32)
                vsSub_64(n, a, b, r);
            else
                vdSub_64(n, a, b, r);
        }
    }
}
