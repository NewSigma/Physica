/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2021-2026 Weibo He
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#pragma once

namespace Physica {
    template<FloatPrec Prec, bool errorTrack>
    Real<Prec> lnGamma(Real<Prec> s) noexcept;

    template<Scalar T>
    T hermiteH(unsigned int n, T x) noexcept {
        if (n == 0)
            return T(1);
        const T double_x = T(2) * x;
        if (n == 1)
            return double_x;

        T old_H = T(1);
        T H = double_x;
        T float_i = T(1);
        for (unsigned int i = 1; i != n; ++i) {
            old_H = fma(double_x, H, T(-2) * float_i * old_H);
            old_H.swap(H);
            float_i += T(1);
        }
        return H;
    }
}
