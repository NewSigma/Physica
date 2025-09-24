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

#include "../Add.h"

namespace Physica {
    template<Vector V1, Vector V2>
    void VectorExpr<ExprType::Add, V1, V2>::assign_mkl(Vector auto& v) const noexcept {
        using Tm = std::conditional<isComplex, typename Tc::MKL_Complex, typename T::MachineType>::type;
        Base::getLHS().assert_assign_mkl(v);
        Base::getRHS().assert_assign_mkl(v);

        size_t n = Base::getLength();
        const auto* a = reinterpret_cast<const Tm*>(Base::getLHS().data());
        const auto* b = reinterpret_cast<const Tm*>(Base::getRHS().data());
        auto* r = reinterpret_cast<Tm*>(v.data());
        if constexpr (isComplex) {
            if constexpr (T::Prec == Float32)
                vcAdd_64(n, a, b, r);
            else
                vzAdd_64(n, a, b, r);
        }
        else {
            if constexpr (T::Prec == Float32)
                vsAdd_64(n, a, b, r);
            else
                vdAdd_64(n, a, b, r);
        }
    }
}
