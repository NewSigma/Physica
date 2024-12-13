/*
 * Copyright 2024 Weibo He.
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

#include <mkl_cblas.h>
#include "../InnerDot.h"

namespace Physica::Core {
    template<Vector T1, Vector T2>
    InnerDot<T1, T2>::ScalarType InnerDot<T1, T2>::calc_mkl() const {
        using MachineType = ScalarType::MachineType;
        const auto* p1 = reinterpret_cast<const MachineType*>(v1.data());
        const auto* p2 = reinterpret_cast<const MachineType*>(v2.data());
        if constexpr (ScalarType::isComplex) {
            ScalarType result;
            if constexpr (ScalarType::Option == Float32)
                cblas_cdotu_sub_64(v1.getLength(), p1, 1, p2, 1, &result);
            else
                cblas_zdotu_sub_64(v1.getLength(), p1, 1, p2, 1, &result);
            return result;
        }
        else {
            if constexpr (ScalarType::Option == Float32)
                return ScalarType(cblas_sdot_64(v1.getLength(), p1, 1, p2, 1));
            else
                return ScalarType(cblas_ddot_64(v1.getLength(), p1, 1, p2, 1));
        }
    }
}
