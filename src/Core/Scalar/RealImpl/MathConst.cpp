/*
 * This file is part of Physica.
 *
 * SPDX-FileCopyrightText: 2024-2026 Weibo He
 * SPDX-License-Identifier: MPL-2.0
 */
#include "Physica/Core/Scalar/Real.h"

using namespace Physica;
/**
 * Consts that need some calculates.
 * Should call new to Const_1 so as to make calculates available.
 */
MathConst<Real<FloatMP>>::MathConst() {
    // 0.31 is the big approximation of ln(2) / ln(10)
    PI = T(calcPI(static_cast<int>(static_cast<double>(MPUnitWidth) * Real<FloatMP>::GlobalPrecision * 0.31) + 1));
    E = T(exp(BasicConst::getInstance()._1));

    PI_2 = T(PI >> 1);
    Minus_PI_2 = T(-PI_2);
}
/**
 * precision is the number of effective digits in decimal.
 * Reference:
 * [1] http://www.pi314.net/eng/salamin.php
 * [2] https://blog.csdn.net/liangbch/article/details/78724041
 */
auto MathConst<Real<FloatMP>>::calcPI(int precision) -> T {
    const auto& basicConst = BasicConst::getInstance();

    T a(static_cast<SignedMPUnit>(1));
    T x(static_cast<SignedMPUnit>(1));
    T b(reciprocal(sqrt(basicConst._2)));
    T c(reciprocal(basicConst._4));

    int goal = 1;
    while (goal < precision) {
        Real y(a);
        a = (a + b) >> 1;
        b = sqrt(b * y);
        y -= a;
        c -= y * y * x;
        x *= basicConst._2;
        goal *= 2;
    }
    a = (a + b) >> 1;
    return a * a / c;
}

auto MathConst<Real<FloatMP>>::getInstance() noexcept -> const This& {
    static MathConst mathConst{};
    return mathConst;
}
