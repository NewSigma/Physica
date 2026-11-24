/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#pragma once

#include <mkl_cblas.h>
#include "Dot.h"

namespace Physica {
    template<Vector LHS, Vector RHS>
    auto Dot<LHS, RHS>::calc_mkl() const noexcept -> T {
        using Tm = decltype(std::declval<T>().toMKL());
        const auto* p1 = reinterpret_cast<const Tm*>(lhs.data_handle());
        const auto* p2 = reinterpret_cast<const Tm*>(rhs.data_handle());
        const auto inc1 = static_cast<MKL_INT64>(lhs.getStride());
        const auto inc2 = static_cast<MKL_INT64>(rhs.getStride());
        if constexpr (T::isComplex()) {
            T result;
            if constexpr (T::Prec == Float32)
                cblas_cdotu_sub_64(lhs.getLength(), p1, inc1, p2, inc2, &result);
            else
                cblas_zdotu_sub_64(lhs.getLength(), p1, inc1, p2, inc2, &result);
            return result;
        }
        else {
            if constexpr (T::Prec == Float32)
                return T(cblas_sdot_64(lhs.getLength(), p1, inc1, p2, inc2));
            else
                return T(cblas_ddot_64(lhs.getLength(), p1, inc1, p2, inc2));
        }
    }
}
